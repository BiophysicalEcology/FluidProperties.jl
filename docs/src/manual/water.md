# Properties of liquid water

The physical properties of liquid water at a given temperature are computed by [`water_properties`](@ref), which returns a
[`WaterProperties`](@ref).

These properties are based on regressions obtained using Grapher from Golden Software on data from Ede, "An Introduction of Heat
Transfer Principles and Calculations," Pergamon Press, 1967, p. 262. The regression was performed by W. Porter on 14 July 1988.
The regressions are valid for temperatures up to 60 °C. Temperatures above 60 °C are clamped to 60 °C for the density.

| Field | Property | Units |
| :---- | :------- | :---- |
| `density` | Density of water | kg m⁻³ |
| `specific_heat` | Specific heat capacity of water | J kg⁻¹ K⁻¹ |
| `thermal_conductivity` | Thermal conductivity of water | W m⁻¹ K⁻¹ |
| `dynamic_viscosity` | Dynamic viscosity of water | kg m⁻¹ s⁻¹ |

```@setup water
using Main.FigureHelpers
using CairoMakie, FluidProperties, Unitful
```

```@example water
using CairoMakie, FluidProperties, Unitful

water_properties(20.0u"°C")
```

Temperatures can be in any unit:

```@example water
water_properties(293.15u"K") == water_properties(20.0u"°C")
```

## Figures

```@example water
twater = collect(0:5:60) .* u"°C" # sequence of water temperatures
water = water_properties.(twater) # compute values

fig = Figure(size = (800, 700))
properties = (
    (:density, u"kg/m^3", "Density (kg m⁻³)"),
    (:specific_heat, u"J/kg/K", "Specific heat (J kg⁻¹ K⁻¹)"),
    (:thermal_conductivity, u"W/m/K", "Thermal conductivity (W m⁻¹ K⁻¹)"),
    (:dynamic_viscosity, u"kg/m/s", "Dynamic viscosity (kg m⁻¹ s⁻¹)"),
)
for (i, (name, unit, label)) in enumerate(properties)
    ax = Axis(fig[fld1(i, 2), mod1(i, 2)]; xlabel = "Water temperature (°C)", ylabel = label)
    lines!(ax, ustrip.(twater), ustrip.(unit, getproperty.(water, name)); linewidth = 2)
end
fig
```
