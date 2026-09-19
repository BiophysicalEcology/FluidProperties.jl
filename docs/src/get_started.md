# Get started

FluidProperties.jl calculates physical properties of air, water vapour and liquid water. Inputs and outputs are
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities.

```julia
using Pkg
Pkg.add("FluidProperties")
```

## Dry air

[`dry_air_properties`](@ref) takes the air temperature and pressure:

```@example get_started
using FluidProperties, Unitful

air = dry_air_properties(20.0u"°C", 101325.0u"Pa")
```

The fields of the result are the properties:

```@example get_started
air.density
```

```@example get_started
uconvert(u"μPa*s", air.dynamic_viscosity)
```

## Humid air

[`wet_air_properties`](@ref) also takes the relative humidity, as a fraction:

```@example get_started
wet_air_properties(20.0u"°C", 0.5, 101325.0u"Pa")
```

## Atmospheric pressure

The pressure at a given elevation can be calculated with [`atmospheric_pressure`](@ref):

```@example get_started
atmospheric_pressure(1000.0u"m")
```

## Vapour pressure

[`vapour_pressure`](@ref) uses the Goff-Gratch equation unless another equation is given:

```@example get_started
vapour_pressure(20.0u"°C")
```

```@example get_started
vapour_pressure(Huang(), 20.0u"°C")
```

## Wet bulb temperature

[`wet_bulb_temperature`](@ref) takes the air temperature, relative humidity and pressure:

```@example get_started
uconvert(u"°C", wet_bulb_temperature(30.0u"°C", 0.5, 101325.0u"Pa"))
```

## Broadcasting

All the functions can be broadcast over arrays, including rasters, see [Application to rasters](tutorials/rasters.md):

```@example get_started
temperatures = collect(0:10:40) .* u"°C"
uconvert.(u"°C", wet_bulb_temperature.(temperatures, 0.5, 101325.0u"Pa"))
```

## Missing values

Missing values are propagated:

```@example get_started
vapour_pressure(missing)
```
