# FluidProperties

[![CI](https://github.com/BiophysicalEcology/FluidProperties.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BiophysicalEcology/FluidProperties.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/BiophysicalEcology/FluidProperties.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BiophysicalEcology/FluidProperties.jl/tree/main)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://biophysicalecology.github.io/FluidProperties.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://biophysicalecology.github.io/FluidProperties.jl/dev/)
[![Docs Build](https://github.com/BiophysicalEcology/FluidProperties.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/BiophysicalEcology/FluidProperties.jl/actions/workflows/Documenter.yml)

Functions to compute properties of air and water. See the [documentation](https://biophysicalecology.github.io/FluidProperties.jl/dev/).

Inputs and outputs are [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities, so any
compatible units can be used. Relative humidity is a fraction between 0 and 1.

```julia
using FluidProperties, Unitful
```

## Dry and humid air

`dry_air_properties` and `wet_air_properties` return the density, viscosity, thermal conductivity, vapour
diffusivity and other properties of air, from the Smithsonian Meteorological Tables.

```julia
dry_air_properties(20.0u"°C", 101325.0u"Pa")
wet_air_properties(20.0u"°C", 0.5, 101325.0u"Pa")
```

`atmospheric_pressure` gives the pressure at an elevation, and `water_properties` the properties of liquid water:

```julia
atmospheric_pressure(1000.0u"m")
water_properties(20.0u"°C")
```

## Vapour pressure

`vapour_pressure` uses the Goff-Gratch equation unless another equation is given: `GoffGratch`, `Teten`, `Huang`, `Bolton`
or a `VapourPressureLookup` table.

```julia
vapour_pressure(20.0u"°C")
vapour_pressure(Huang(), 20.0u"°C")
```

## Wet bulb temperature

`wet_bulb_temperature` takes air temperature, fractional relative humidity and pressure, and a
`WetBulbMethod`. The default is `DaviesJones()`.

```julia
wet_bulb_temperature(30.0u"°C", 0.5, 101325.0u"Pa")
wet_bulb_temperature(Stull(), 30.0u"°C", 0.5, 101325.0u"Pa")
wet_bulb_temperature(Smithsonian(; vapour_pressure_equation=Huang()), 30.0u"°C", 0.5, 101325.0u"Pa")
```

| Method | Basis |
|---|---|
| `DaviesJones` | Pseudoadiabatic wet bulb (Davies-Jones 2008); accepts `SpecificHumidity` |
| `Stull` | Empirical fit (Stull 2011); ignores pressure |
| `Barenbrug` | Psychrometric equation (Barenbrug 1947) |
| `Smithsonian` | Psychrometric equation (List 1971) |
| `EnergyBalance` | Adiabatic saturation energy balance (Zhang et al. 2021) |

`wet_bulb_properties` also returns the equivalent and equivalent potential temperatures for
`DaviesJones`. `natural_wet_bulb_temperature` adjusts a psychrometric wet bulb for wind speed.

## Broadcasting

All the functions broadcast over arrays, including rasters:

```julia
temperatures = collect(0:10:40) .* u"°C"
wet_bulb_temperature.(temperatures, 0.5, 101325.0u"Pa")
```

See the [documentation](https://biophysicalecology.github.io/FluidProperties.jl/dev/) for the manual, figures of the
properties of air and an example applying the functions to rasters.
