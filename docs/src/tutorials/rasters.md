# Application to rasters

The functions of FluidProperties.jl can be applied to gridded datasets (rasters) by broadcasting, so that maps of
the properties of air can be made with [Rasters.jl](https://github.com/rafaqz/Rasters.jl). The data are read with
[RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl).

The example below adjusts relative humidity for a change in temperature, assuming the absolute humidity (vapor
density) remains constant, and then maps the density of the air and the wet bulb temperature.

## Setup

Install the packages used in this tutorial:

```julia
using Pkg
Pkg.add(["FluidProperties", "Unitful", "Rasters", "RasterDataSources", "NCDatasets", "CairoMakie"])
```

To download data you will need to specify a folder to put it in. You can do this by assigning the environment variable RASTERDATASOURCES_PATH:

```julia
ENV["RASTERDATASOURCES_PATH"] = joinpath(homedir(), "RasterDataSources") # or "/your/path/here"
```

## Acquiring the data

We'll use the CRU CL 2.0 dataset (New et al. 2002), monthly mean climate for the global land surface at 10 minute
resolution for 1961-1990, which is available through RasterDataSources.jl as `CRUCL2`. It contains the mean temperature
(`tmp`), the diurnal temperature range (`dtr`), the relative humidity (`reh`), wind speed (`wnd`) and elevation (`elv`),
with a month dimension `Ti` for all but elevation.

```@example rasters
using Rasters, RasterDataSources, NCDatasets
using FluidProperties, Unitful
using CairoMakie

path = getraster(CRUCL2)
climate = RasterStack(path; lazy = true)
```

The data are opened lazily, so nothing is read into memory yet. Select January, and read it.

```@example rasters
january = read(climate[Ti(1)])
```

The maximum and minimum air temperatures are the mean temperature plus and minus half of the diurnal range,
and the relative humidity is a percentage, which we convert to a fraction. We give the rasters units, so that all
the calculations are checked for unit mistakes.

```@example rasters
Tmean = january.tmp .* u"°C"
Tmax = (january.tmp .+ january.dtr ./ 2) .* u"°C"
Tmin = (january.tmp .- january.dtr ./ 2) .* u"°C"
RHmean = january.reh ./ 100
elevation = january.elv .* u"m"
nothing # hide
```

Let's plot the January mean air temperature and relative humidity:

```@example rasters
f = Figure(size = (700, 750))
a1 = Axis(f[1, 1]; title = "Mean air temperature, January (°C)")
p1 = heatmap!(a1, ustrip.(Tmean))
Colorbar(f[1, 2], p1)
a2 = Axis(f[2, 1]; title = "Mean relative humidity, January")
p2 = heatmap!(a2, RHmean)
Colorbar(f[2, 2], p2)
f
```

## Relative humidity at the minimum temperature

In the calculations below, the vapor pressure is calculated from the mean relative humidity at the mean air temperature
with [`vapour_pressure`](@ref) (broadcasting over the raster), then the saturation vapor pressure at the minimum air
temperature is obtained and the relative humidity at the minimum temperature is computed from the ratio of actual to saturation vapor pressure.

```@example rasters
e = RHmean .* vapour_pressure.(Tmean) # vapor pressure
esat = vapour_pressure.(Tmin) # saturation vapor pressure at the minimum temperature

# compute new relative humidity for minimum air temperature
RHmax = e ./ esat
# conditional replace of any values >1 with 1
RHmax = min.(RHmax, 1.0)
nothing # hide
```

Now plot the results:

```@example rasters
f = Figure(size = (700, 1100))
a1 = Axis(f[1, 1]; title = "Vapor pressure, January (Pa)")
p1 = heatmap!(a1, ustrip.(u"Pa", e))
Colorbar(f[1, 2], p1)
a2 = Axis(f[2, 1]; title = "Saturation vapor pressure at Tmin (Pa)")
p2 = heatmap!(a2, ustrip.(u"Pa", esat))
Colorbar(f[2, 2], p2)
a3 = Axis(f[3, 1]; title = "Relative humidity at Tmin, January")
p3 = heatmap!(a3, RHmax; colorrange = (0, 1))
Colorbar(f[3, 2], p3)
f
```

## Pressure and density of the air

The atmospheric pressure is calculated from the elevation with [`atmospheric_pressure`](@ref), and the properties of
humid air at the mean temperature and humidity with [`wet_air_properties`](@ref). A `WetAirProperties` is returned for every
cell, and the property we want is taken from each with `map`, leaving missing cells missing.

```@example rasters
pressure = atmospheric_pressure.(elevation)
air = wet_air_properties.(Tmean, RHmean, pressure)
density = map(a -> ismissing(a) ? missing : a.density, air)

f = Figure(size = (700, 750))
a1 = Axis(f[1, 1]; title = "Atmospheric pressure (Pa)")
p1 = heatmap!(a1, ustrip.(u"Pa", pressure))
Colorbar(f[1, 2], p1)
a2 = Axis(f[2, 1]; title = "Density of humid air, January (kg m⁻³)")
p2 = heatmap!(a2, ustrip.(u"kg/m^3", density))
Colorbar(f[2, 2], p2)
f
```

## Wet bulb temperature

[`wet_bulb_temperature`](@ref) also broadcasts over rasters. Here it is calculated for the January mean conditions
with the default [`DaviesJones`](@ref) method, which depends on the pressure, and with the empirical [`Stull`](@ref) formula.

```@example rasters
wetbulb = wet_bulb_temperature.(Tmean, RHmean, pressure)
wetbulb_stull = wet_bulb_temperature.(Stull(), Tmean, RHmean, pressure)
nothing # hide
```

```@example rasters
f = Figure(size = (700, 750))
a1 = Axis(f[1, 1]; title = "Wet bulb temperature, Davies-Jones (°C)")
p1 = heatmap!(a1, ustrip.(u"°C", wetbulb))
Colorbar(f[1, 2], p1)
a2 = Axis(f[2, 1]; title = "Stull minus Davies-Jones (K)")
p2 = heatmap!(a2, ustrip.(u"K", wetbulb_stull .- wetbulb); colorrange = (-1, 1), colormap = :balance)
Colorbar(f[2, 2], p2)
f
```

## Timing

The calculations above are broadcasts over every land cell of a global raster at 10 minute resolution. The time taken
by each function, after compilation, gives an idea of the speed. The functions are type-stable and do not allocate, so
the time is proportional to the number of cells.

```@example rasters
using Markdown

ncells = count(!ismissing, Tmean)
timings = (
    atmospheric_pressure = @elapsed(atmospheric_pressure.(elevation)),
    vapour_pressure = @elapsed(vapour_pressure.(Tmean)),
    wet_air_properties = @elapsed(wet_air_properties.(Tmean, RHmean, pressure)),
    wet_bulb_temperature_DaviesJones = @elapsed(wet_bulb_temperature.(Tmean, RHmean, pressure)),
    wet_bulb_temperature_Stull = @elapsed(wet_bulb_temperature.(Stull(), Tmean, RHmean, pressure)),
)

rows = ["| Function | Time (s) | Time per cell (ns) |", "| :------- | -------: | -----------------: |"]
for (name, seconds) in pairs(timings)
    push!(rows, "| `$name` | $(round(seconds; sigdigits = 2)) | $(round(Int, 1e9 * seconds / ncells)) |")
end
Markdown.parse("$ncells cells:\n\n" * join(rows, "\n"))
```

## References

New M, Lister D, Hulme M, Makin I. 2002. A high-resolution data set of surface climate over global land areas. Climate Research 21: 1-25.
