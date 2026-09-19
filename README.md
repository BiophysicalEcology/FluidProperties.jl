# FluidProperties

[![CI](https://github.com/BiophysicalEcology/FluidProperties.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BiophysicalEcology/FluidProperties.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/BiophysicalEcology/FluidProperties.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BiophysicalEcology/FluidProperties.jl/tree/main)

Functions to compute properties of air and water

## Wet bulb temperature

`wet_bulb_temperature` takes air temperature, fractional relative humidity and pressure, and a
`WetBulbMethod`. The default is `DaviesJones()`.

```julia
using FluidProperties, Unitful
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
| `EnergyBalance` | Adiabatic saturation energy balance |

`wet_bulb_properties` also returns the equivalent and equivalent potential temperatures for
`DaviesJones`. `natural_wet_bulb_temperature` adjusts a psychrometric wet bulb for wind speed.
