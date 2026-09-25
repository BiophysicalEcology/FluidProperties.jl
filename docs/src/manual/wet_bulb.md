# Wet bulb temperature

The wet bulb temperature is the temperature read from a thermometer with its bulb covered with a wet wick
and ventilated. It is a function of air temperature, humidity and pressure, and is calculated by
[`wet_bulb_temperature`](@ref). The relationship of wet bulb temperature to dry bulb temperature, vapor density,
and relative humidity is shown in the psychometric diagram in [Figure 15](humid_air.md#Figure-15.-Vapour-density-over-water).

```@setup wetbulb
using Main.FigureHelpers
using CairoMakie, FluidProperties, Unitful
```

```@example wetbulb
using CairoMakie, FluidProperties, Unitful

wet_bulb_temperature(30.0u"°C", 0.5, 101325.0u"Pa")
```

The wet bulb temperature is returned in kelvin. Air temperature and pressure can be given in any unit, and
relative humidity is a fraction between 0 and 1. A humidity outside of this range throws a `DomainError`.

```@example wetbulb
uconvert(u"°C", wet_bulb_temperature(86.0u"°F", 0.5, 1013.25u"hPa"))
```

## Methods

The formulation is chosen by passing a [`WetBulbMethod`](@ref) as the first argument. The default is [`DaviesJones`](@ref).

| Method | Basis |
| :----- | :---- |
| [`DaviesJones`](@ref) | Pseudoadiabatic wet bulb, Davies-Jones (2008); accepts [`SpecificHumidity`](@ref) |
| [`Stull`](@ref) | Empirical fit, Stull (2011); ignores pressure |
| [`Barenbrug`](@ref) | Psychrometric equation, Barenbrug (1947) |
| [`Smithsonian`](@ref) | Psychrometric equation, List (1971) |
| [`EnergyBalance`](@ref) | Adiabatic saturation energy balance, Zhang et al. (2021) |

```@example wetbulb
T = 30.0u"°C"
RH = 0.5
P = 101325.0u"Pa"

methods = (DaviesJones(), Stull(), Barenbrug(), Smithsonian(), EnergyBalance())
[nameof(typeof(m)) => uconvert(u"°C", wet_bulb_temperature(m, T, RH, P)) for m in methods]
```

The psychrometric methods, [`Barenbrug`](@ref), [`Smithsonian`](@ref) and [`EnergyBalance`](@ref), solve an equation for the wet bulb
temperature and take a `vapour_pressure_equation` ([`GoffGratch`](@ref) by default) and the solver settings `tolerance` and `max_iterations`.
With the [`Bolton`](@ref) vapour pressure equation, [`Smithsonian`](@ref) is the equation used in the NOAA wet bulb calculator
(Lemke and Kjellstrom 2012).

```@example wetbulb
wet_bulb_temperature(Smithsonian(; vapour_pressure_equation = Bolton()), T, RH, P)
```

### Comparison of the methods

```@example wetbulb
tairs = collect(0:2:50) .* u"°C" # sequence of air temperatures
twb(method, rh = 0.5, p = 101325.0u"Pa") = uconvert.(u"°C", wet_bulb_temperature.(method, tairs, rh, p))

fig, ax = figure_axis("Air temperature (°C)", "Wet bulb temperature (°C)")
for method in methods
    lines!(ax, ustrip.(tairs), ustrip.(twb(method)); linewidth = 2, label = string(nameof(typeof(method))))
end
axislegend(ax; position = :lt)
fig
```

The methods differ from [`DaviesJones`](@ref) by less than about 1 K at 50 °C and 50 % relative humidity, and
by less at lower temperatures.

```@example wetbulb
fig, ax = figure_axis("Air temperature (°C)", "Difference from Davies-Jones (K)")
reference = twb(DaviesJones())
for method in methods[2:end]
    lines!(ax, ustrip.(tairs), ustrip.(u"K", twb(method) .- reference); linewidth = 2, label = string(nameof(typeof(method))))
end
axislegend(ax; position = :lb)
fig
```

### Relative humidity and pressure

The wet bulb temperature rises with relative humidity and is the air temperature at saturation. At lower pressure
it is lower, as there is less air to be cooled by evaporation.

```@example wetbulb
fig, ax = figure_axis("Air temperature (°C)", "Wet bulb temperature (°C)")
for rh in (0.1, 0.25, 0.5, 0.75, 1.0)
    lines!(ax, ustrip.(tairs), ustrip.(twb(DaviesJones(), rh)); linewidth = 2, label = "$(round(Int, 100rh)) % rh")
end
axislegend(ax; position = :lt)
fig
```

```@example wetbulb
fig, ax = figure_axis("Air temperature (°C)", "Wet bulb temperature (°C)")
for p in (101325.0u"Pa", 85000.0u"Pa", 70000.0u"Pa")
    lines!(ax, ustrip.(tairs), ustrip.(twb(DaviesJones(), 0.5, p)); linewidth = 2, label = "$p")
end
axislegend(ax; position = :lt)
fig
```

## Davies-Jones method

[`wet_bulb_properties`](@ref) returns the equivalent temperature and equivalent potential temperature as well as the
wet bulb temperature, as a [`WetBulbProperties`](@ref).
The lifting condensation temperature (Bolton 1980, eqn 22), the moist potential temperature (eqn 24) and the equivalent
potential temperature (eqn 39) are found first, then the wet bulb temperature by Newton-Raphson iteration from the
first guess of Davies-Jones (2008). The wet bulb and equivalent temperatures are `NaN` where the equivalent temperature
is outside 200-600 K.

```@example wetbulb
wet_bulb_properties(30.0u"°C", 0.5, 101325.0u"Pa")
```

The humidity can also be the specific humidity, which is treated as a mixing ratio:

```@example wetbulb
wet_bulb_properties(30.0u"°C", SpecificHumidity(0.0135), 101325.0u"Pa").wet_bulb_temperature
```

The iteration is chosen with the `convergence` keyword: [`RefinedNewton`](@ref) (the default) iterates to a tolerance
and [`FixedNewton`](@ref) does at most four iterations.

```@example wetbulb
wet_bulb_temperature(DaviesJones(; convergence = FixedNewton()), T, RH, P) - wet_bulb_temperature(T, RH, P)
```

## Natural wet bulb temperature

The natural wet bulb temperature is the temperature of a wetted thermometer bulb in the natural
environment of wind and sun. [`natural_wet_bulb_temperature`](@ref) calculates it from the psychrometric wet bulb temperature
and the wind speed (Bernard and Pourmoghani 1999). If the globe temperature is more than 4 K above the air temperature
the effect of radiant heat is included. The air temperature is the default globe temperature.

```@example wetbulb
tw = wet_bulb_temperature(30.0u"°C", 0.5, 101325.0u"Pa")
natural_wet_bulb_temperature(30.0u"°C", tw, 1.0u"m/s")
```

```@example wetbulb
winds = (0.01:0.01:4.0) .* u"m/s" # sequence of wind speeds

fig, ax = figure_axis("Wind speed (m s⁻¹)", "Natural wet bulb temperature (°C)")
lines!(ax, ustrip.(winds), ustrip.(uconvert.(u"°C", natural_wet_bulb_temperature.(30.0u"°C", tw, winds))); linewidth = 2, label = "no radiant heat")
lines!(ax, ustrip.(winds), ustrip.(uconvert.(u"°C", natural_wet_bulb_temperature.(30.0u"°C", tw, winds; globe_temperature = 40.0u"°C"))); linewidth = 2, label = "globe 10 K above air")
hlines!(ax, ustrip(uconvert(u"°C", tw)); linestyle = :dash, color = :grey, label = "psychrometric")
axislegend(ax; position = :rt)
fig
```

## References

Barenbrug AWT. 1947. Psychrometry and psychrometric charts. Journal of the Chemical, Metallurgical and Mining Society of South Africa, May: 393-417.

Bernard TE, Pourmoghani M. 1999. Prediction of workplace wet bulb global temperature. Applied Occupational and Environmental Hygiene 14: 126-134.

Bolton D. 1980. The computation of equivalent potential temperature. Monthly Weather Review 108: 1046-1053.

Davies-Jones R. 2008. An efficient and accurate method for computing the wet-bulb temperature along pseudoadiabats. Monthly Weather Review 136: 2764-2785.

Lemke B, Kjellstrom T. 2012. Calculating workplace WBGT from meteorological data: a tool for climate change assessment. Industrial Health 50: 267-278.

List RJ. 1971. Smithsonian Meteorological Tables. Smithsonian Institution Press.

Stull R. 2011. Wet-bulb temperature from relative humidity and air temperature. Journal of Applied Meteorology and Climatology 50: 2267-2269.

Zhang Y, Held I, Fueglistaler S. 2021. Projections of tropical heat stress constrained by atmospheric dynamics. Nature Geoscience 14: 133-137.
