# Properties of humid air

The properties of humid air in this section are computed by [`wet_air_properties`](@ref), which takes
air temperature, relative humidity (as a fraction) and pressure, and returns a [`WetAirProperties`](@ref).
Saturation vapour pressure is calculated by [`vapour_pressure`](@ref).

```@setup humidair
using Main.FigureHelpers
using CairoMakie, FluidProperties, Unitful
```

```@example humidair
using CairoMakie, FluidProperties, Unitful

tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
air = wet_air_properties.(tairs, 1.0, 100000.0u"Pa"); # a `WetAirProperties` for each temperature at 100 % rh
air[1]
```

## Figure 3. Dew point temperature

The saturation vapor pressure at the dew point is calculated with the Goff-Gratch equations, which
switch to saturation over ice below freezing. The dashed lines are the Tetens equation
([`Teten`](@ref)) over water, with $\alpha = 7.5$ and $\beta = 237.3$, and

```math
e_d^* = 100 \times 10^{\left[0.7858 + \frac{t_d \alpha}{t_d + \beta}\right]}
```

over ice, with $\alpha = 9.5$ and $\beta = 265.5$.

```@example humidair
dps_water = collect(0:5:50) .* u"°C" # sequence of dew points over water
dps_ice = collect(-50:5:0) .* u"°C" # sequence of dew points over ice
esat_water = vapour_pressure.(dps_water) # compute values
esat_ice = vapour_pressure.(dps_ice)

# Tetens equation
tetens_water = vapour_pressure.(Teten(), dps_water)
tetens_ice = 100 .* 10 .^ (0.7858 .+ ustrip.(dps_ice) .* 9.5 ./ (ustrip.(dps_ice) .+ 265.5)) .* u"Pa"

fig = Figure(size = (700, 800))
ax1 = Axis(fig[1, 1]; xlabel = "Saturation vapor pressure at tᵈ over water (Pa)", ylabel = "Dew point temperature (°C)")
lines!(ax1, ustrip.(u"Pa", esat_water), ustrip.(dps_water); linewidth = 2, label = "Goff-Gratch")
lines!(ax1, ustrip.(u"Pa", tetens_water), ustrip.(dps_water); linewidth = 2, linestyle = :dash, label = "Tetens")
axislegend(ax1; position = :rb)
ax2 = Axis(fig[2, 1]; xlabel = "Saturation vapor pressure at tᵈ over ice (Pa)", ylabel = "Dew point temperature (°C)")
lines!(ax2, ustrip.(u"Pa", esat_ice), ustrip.(dps_ice); linewidth = 2, label = "Goff-Gratch")
lines!(ax2, ustrip.(u"Pa", tetens_ice), ustrip.(dps_ice); linewidth = 2, linestyle = :dash, label = "Tetens")
axislegend(ax2; position = :rb)
fig
```

## Figure 9. Mixing ratio over water at 100 000 Pa

```math
r_w = \frac{0.62570 \, e}{p - 1.0060 \, e}
```

```@example humidair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
rw(rh) = getproperty.(wet_air_properties.(tairs, rh, 100000.0u"Pa"), :mixing_ratio) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Mixing ratio over water at 100 000 Pa (kg kg⁻¹)")
for rh in (1.0, 0.5, 0.25)
    lines!(ax, ustrip.(tairs), rw(rh); linewidth = 2, label = "$(round(Int, 100rh)) % rh")
end
axislegend(ax; position = :lt)
fig
```

## Figure 10. Specific heat of air at 100 000 Pa

```math
c_p = \frac{1004.84 + (1864.40 \, r_w)}{1 + r_w}
```

```@example humidair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
cp(rh) = getproperty.(wet_air_properties.(tairs, rh, 100000.0u"Pa"), :specific_heat) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Specific heat of air at 100 000 Pa (J kg⁻¹ K⁻¹)")
for rh in (1.0, 0.5, 0.25)
    lines!(ax, ustrip.(tairs), ustrip.(u"J/kg/K", cp(rh)); linewidth = 2, label = "$(round(Int, 100rh)) % rh")
end
axislegend(ax; position = :lt)
fig
```

## Figure 15. Vapour density over water

```math
\rho_v = \frac{e}{461.5 \, (t + 273.15)}
```

The curves are lines of constant relative humidity. The straight lines are lines of constant wet bulb temperature
(see [Wet bulb temperature](wet_bulb.md)), from the saturated vapour density at the wet bulb temperature to zero vapour
density, found by rearranging the psychrometer equation of the Smithsonian tables ([`Smithsonian`](@ref)).

```@example humidair
tairs = collect(0:5:60) .* u"°C" # sequence of air temperatures at which to obtain vapor density
P = 101325.0u"Pa"
vd(rh) = getproperty.(wet_air_properties.(tairs, rh, P), :vapour_density) # obtain vapor density at each relative humidity

fig, ax = figure_axis("Air temperature (°C)", "Vapor density over water (kg m⁻³)"; limits = (nothing, (0, 0.09)))
for rh in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
    lines!(ax, ustrip.(tairs), ustrip.(u"kg/m^3", vd(rh)); linewidth = 2, color = :black)
end

# get dry bulb temperatures for a set of saturated wet bulb temperatures
for wb in collect(0:2:50) .* u"°C" # sequence of wet bulb temperatures
    esat = vapour_pressure(wb)
    # rearrangement of the psychrometer equation
    db = wb + esat / (0.00066 / u"K" * (1 + 0.00115 / u"K" * (uconvert(u"K", wb) - 273.15u"K")) * P)
    vd_wb = wet_air_properties(wb, 1.0, P).vapour_density
    lines!(ax, [ustrip(wb), ustrip(db)], [ustrip(u"kg/m^3", vd_wb), 0.0]; linewidth = 1, color = :steelblue)
end
fig
```

## Figure 16. Vapour pressure over water

```math
e = 461.50 \, \rho_v (t + 273.15)
```

```@example humidair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
e(rh) = getproperty.(wet_air_properties.(tairs, rh, 101325.0u"Pa"), :vapour_pressure) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Vapor pressure over water (Pa)"; limits = (nothing, (0, 10e3)))
for rh in (1.0, 0.5, 0.25)
    lines!(ax, ustrip.(tairs), ustrip.(u"Pa", e(rh)); linewidth = 2, label = "$(round(Int, 100rh)) % rh")
end
axislegend(ax; position = :lt)
fig
```

## Figure 17. Virtual temperature increment at 100 000 Pa

```math
\Delta T_v = T \left[\frac{1 + \frac{r_w}{0.622}}{1 + r_w}\right] - T, \qquad T = t + 273.15
```

```@example humidair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
tvinc(rh) = getproperty.(wet_air_properties.(tairs, rh, 100000.0u"Pa"), :virtual_temp_increment) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Virtual temperature increment at 100 000 Pa (K)")
for rh in (1.0, 0.5, 0.25)
    lines!(ax, ustrip.(tairs), ustrip.(u"K", tvinc(rh)); linewidth = 2, label = "$(round(Int, 100rh)) % rh")
end
axislegend(ax; position = :lt)
fig
```

## Figure 18. Water potential

```math
\psi = (4.615 \times 10^{5}) (t + 273.15) \ln\left(\frac{rh}{100}\right)
```

```@example humidair
rhs = 0.65:0.05:1.0 # relative humidities (fractions)
wtrpot(t) = getproperty.(wet_air_properties.(t, rhs, 101325.0u"Pa"), :water_potential) # compute values

fig, ax = figure_axis("Relative humidity (%)", "-1 × water potential (Pa)")
for t in collect(0:10:50) .* u"°C"
    lines!(ax, 100 .* rhs, -ustrip.(u"Pa", wtrpot(t)); linewidth = 2, label = "$t")
end
axislegend(ax; position = :lt)
fig
```

Near saturation the water potential is small:

```@example humidair
rhs = 0.993:0.0005:1.0 # relative humidities (fractions)

fig, ax = figure_axis("Relative humidity (%)", "-1 × water potential (Pa)")
for t in collect(0:10:50) .* u"°C"
    lines!(ax, 100 .* rhs, -ustrip.(u"Pa", wtrpot(t)); linewidth = 2, label = "$t")
end
axislegend(ax; position = :lt)
fig
```
