# Saturation vapour pressure

Saturation vapour pressure is calculated by [`vapour_pressure`](@ref) for a given air temperature. All humidity
calculations in [`wet_air_properties`](@ref) are made with the internationally accepted Goff-Gratch equation
([`GoffGratch`](@ref)), which is the default. A mathematically less cumbersome and often more useful alternative
was reported by Tetens (Murray 1976), available as [`Teten`](@ref); calculations with this equation differ
negligibly from those with the Goff-Gratch expression.

```@setup vapourpressure
using Main.FigureHelpers
using CairoMakie, FluidProperties, Unitful
```

| Equation | Description |
| :------- | :---------- |
| [`GoffGratch`](@ref) | Goff-Gratch equations over ice below 0 °C and water above. The default. |
| [`Teten`](@ref) | Tetens equation. Low accuracy but very fast, with a single `exp` call. |
| [`Huang`](@ref) | Huang (2018). High accuracy from -100 to 100 °C and reasonable performance. |
| [`Bolton`](@ref) | Bolton (1980) eqn 10, over liquid water. Used by [`DaviesJones`](@ref) for the wet bulb. |
| [`ClausiusClapeyron`](@ref) | Integrated Clausius-Clapeyron relation, over ice below the triple point and water above. Matches [Thermodynamics.jl](https://github.com/CliMA/Thermodynamics.jl). |
| [`VapourPressureLookup`](@ref) | A lookup table with linear interpolation, with no transcendental function calls. |

```@example vapourpressure
using CairoMakie, FluidProperties, Unitful

vapour_pressure(20.0u"°C") # Goff-Gratch by default
```

```@example vapourpressure
vapour_pressure(Huang(), 20.0u"°C")
```

Temperatures may be given in any unit:

```@example vapourpressure
vapour_pressure(Teten(), 293.15u"K")
```

## Comparison of the equations

Each equation is compared with the Goff-Gratch equation. The lookup table is built from Goff-Gratch with a step of 0.5 K here.

```@example vapourpressure
tairs = collect(-20:1:50) .* u"°C" # sequence of air temperatures
lookup = VapourPressureLookup(GoffGratch(); tmin = -40.0u"°C", tmax = 60.0u"°C", step = 0.5u"K")
equations = (GoffGratch = GoffGratch(), Teten = Teten(), Huang = Huang(), Bolton = Bolton(),
    ClausiusClapeyron = ClausiusClapeyron(), Lookup = lookup)

fig, ax = figure_axis("Air temperature (°C)", "Saturation vapor pressure (Pa)")
for (name, equation) in pairs(equations)
    lines!(ax, ustrip.(tairs), ustrip.(u"Pa", vapour_pressure.(equation, tairs)); linewidth = 2, label = string(name))
end
axislegend(ax; position = :lt)
fig
```

```@example vapourpressure
reference = vapour_pressure.(GoffGratch(), tairs)

fig, ax = figure_axis("Air temperature (°C)", "Difference from Goff-Gratch (%)")
for (name, equation) in pairs(equations)
    name == :GoffGratch && continue
    difference = 100 .* (vapour_pressure.(equation, tairs) ./ reference .- 1)
    lines!(ax, ustrip.(tairs), difference; linewidth = 2, label = string(name))
end
axislegend(ax; position = :lt)
fig
```

Bolton's equation is for liquid water, so it departs from the Goff-Gratch equation over ice below 0 °C.

## Clausius-Clapeyron and Thermodynamics.jl

The other equations are empirical fits. [`ClausiusClapeyron`](@ref) instead integrates the Clausius-Clapeyron relation
from the triple point of water, assuming the specific heats of vapour, liquid water and ice are constant (the
Rankine-Kirchhoff approximation), so that the latent heat varies linearly with temperature:

```math
e^* = p_{tr} \left(\frac{T}{T_{tr}}\right)^{Δc_p/R_v}
      \exp\left[\frac{L_0 - Δc_p T_0}{R_v}\left(\frac{1}{T_{tr}} - \frac{1}{T}\right)\right]
```

This is the formulation used by [Thermodynamics.jl](https://github.com/CliMA/Thermodynamics.jl), the moist
thermodynamics package of the Climate Modeling Alliance (CliMA). The default constants are those of
[ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl), so the results match Thermodynamics.jl over liquid water
above the triple point and over ice below it. The constants are keywords, and `ice = false` gives saturation over
supercooled water at all temperatures:

```@example vapourpressure
vapour_pressure(ClausiusClapeyron(; ice = false), -10.0u"°C")
```

With its default parameters, Thermodynamics.jl's `saturation_vapor_pressure(param_set, T)` blends the liquid and ice
curves between -40 °C and 0 °C with a temperature-dependent liquid fraction (Pressel et al. 2015), so it differs from
`ClausiusClapeyron()` in that range.

The two packages overlap only in the humidity calculations, and serve different purposes:

| | FluidProperties.jl | Thermodynamics.jl |
| :-- | :-- | :-- |
| Purpose | Properties of air and water for biophysical ecology and microclimate models | Moist thermodynamics for climate and weather models |
| Scope | Transport properties (viscosity, conductivity, diffusivity), liquid water, wet bulb temperatures | Energies, entropies, potential temperatures, phase equilibrium (saturation adjustment) |
| Vapour pressure | Several empirical equations, and the Clausius-Clapeyron relation | Clausius-Clapeyron relation, with liquid-ice mixtures |
| Constants | Those of each source equation | A single parameter set from ClimaParams.jl |
| Inputs | Unitful.jl quantities, `missing` propagated | Plain floats in SI units, GPU and automatic differentiation support |

## Adding an equation

Each equation is a type, and [`vapour_pressure`](@ref) has a method for each. Julia chooses the method from the
types of all the arguments (multiple dispatch), so an equation can be added without changing the package: define a
subtype of [`VapourPressureEquation`](@ref) and a method for it. The Magnus form fitted by Alduchov and Eskridge
(1996), over liquid water:

```@example vapourpressure
struct AlduchovEskridge <: VapourPressureEquation end

function FluidProperties.vapour_pressure(::AlduchovEskridge, T)
    Tc = ustrip(u"°C", T)
    return 6.1094 * exp(17.625 * Tc / (243.04 + Tc)) * 100u"Pa"
end

vapour_pressure(AlduchovEskridge(), 20.0u"°C")
```

It now works wherever the others do, including in [`wet_air_properties`](@ref) and in the packages built on
FluidProperties.jl, such as Microclimate.jl:

```@example vapourpressure
wet_air_properties(20.0u"°C", 0.5, 101325.0u"Pa"; vapour_pressure_equation = AlduchovEskridge()).vapour_pressure
```

## References

Alduchov OA, Eskridge RE. 1996. Improved Magnus form approximation of saturation vapor pressure. Journal of Applied Meteorology 35: 601-609.

Bolton D. 1980. The computation of equivalent potential temperature. Monthly Weather Review 108: 1046-1053.

Huang J. 2018. A simple accurate formula for calculating saturation vapor pressure of water and ice. Journal of Applied Meteorology and Climatology 57: 1265-1272.

Murray BR. 1976. On the computation of saturation vapor pressure. Journal of Applied Meteorology 6:203-204.

Pressel KG, Kaul CM, Schneider T, Tan Z, Mishra S. 2015. Large-eddy simulation in an anelastic framework with closed water and entropy balances. Journal of Advances in Modeling Earth Systems 7: 1425-1456.
