# Properties of dry air

The properties of dry air in this section are computed by [`dry_air_properties`](@ref), which takes
air temperature and pressure and returns a [`DryAirProperties`](@ref). Standard atmospheric pressure is
computed by [`atmospheric_pressure`](@ref) and the latent heat of vaporization of water by
[`enthalpy_of_vaporisation`](@ref).

```@setup dryair
using Main.FigureHelpers
using CairoMakie, FluidProperties, Unitful
```

```@example dryair
using CairoMakie, FluidProperties, Unitful

tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
air = dry_air_properties.(tairs, 101325.0u"Pa"); # a `DryAirProperties` for each temperature
air[1]
```

The properties of the `air` at each temperature are found by broadcasting over the fields, e.g. `getproperty.(air, :density)`.

## Figure 1. Black-body emittance

```math
\phi = 5.67032 \times 10^{-8} (t + 273.15)^{4}
```

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
bbemit = getproperty.(dry_air_properties.(tairs, 101325.0u"Pa"), :blackbody_emission) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Black-body emittance (W m⁻²)"; limits = (nothing, (0, 1000)))
lines!(ax, ustrip.(tairs), ustrip.(u"W/m^2", bbemit); linewidth = 2)
fig
```

## Figure 2. Density of dry air

```math
\rho = \frac{P}{287.04 (t + 273.15)}
```

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
density(P) = getproperty.(dry_air_properties.(tairs, P), :density) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Density of dry air (kg m⁻³)"; limits = (nothing, (1, 1.5)))
for P in (85000.0u"Pa", 100000.0u"Pa")
    lines!(ax, ustrip.(tairs), ustrip.(u"kg/m^3", density(P)); linewidth = 2, label = "$P")
end
axislegend(ax; position = :rt)
fig
```

## Figure 4. Diffusivity of water vapor in air

```math
D = D_0 \left(\frac{T}{T_0}\right)^{n} \frac{\rho}{\rho_0}, \quad D_0 = 2.26 \times 10^{-5}, \quad T = t + 273.15, \quad n = 1.81, \quad \rho_0 = 10^5
```

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
diffusivity(P) = getproperty.(dry_air_properties.(tairs, P), :vapour_diffusivity) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Diffusivity of water vapor in air (m² s⁻¹)")
for P in (70000.0u"Pa", 85000.0u"Pa", 100000.0u"Pa")
    lines!(ax, ustrip.(tairs), ustrip.(u"m^2/s", diffusivity(P)); linewidth = 2, label = "$P")
end
axislegend(ax; position = :lt)
fig
```

## Figure 5. Dynamic viscosity of air

```math
\mu = \mu_0 \left[\frac{T_0 + C}{T + C} \left(\frac{T}{T_0}\right)^{1.5}\right], \quad \mu_0 = 1.8325 \times 10^{-5}, \quad T_0 = 296.16, \quad C = 120, \quad T = t + 273.15
```

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
visdyn = getproperty.(dry_air_properties.(tairs, 101325.0u"Pa"), :dynamic_viscosity) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Dynamic viscosity of air (kg m⁻¹ s⁻¹)")
lines!(ax, ustrip.(tairs), ustrip.(u"kg/m/s", visdyn); linewidth = 2)
fig
```

## Figure 6. Group of variables in Grashof number

```math
\gamma = \frac{g \beta}{\nu^{2}}
```

The Grashof number is obtained by multiplying by $\Delta T L^3$.

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
ggroup = getproperty.(dry_air_properties.(tairs, 101325.0u"Pa"), :grashof_coefficient) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Group of variables in Grashof number (m⁻³ K⁻¹)")
lines!(ax, ustrip.(tairs), ustrip.(u"m^-3/K", ggroup); linewidth = 2)
fig
```

## Figure 7. Kinematic viscosity of air at 100 000 Pa

```math
\nu = \frac{\mu}{\rho}
```

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
viskin = getproperty.(dry_air_properties.(tairs, 100000.0u"Pa"), :kinematic_viscosity) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Kinematic viscosity of air at 100 000 Pa (m² s⁻¹)")
lines!(ax, ustrip.(tairs), ustrip.(u"m^2/s", viskin); linewidth = 2)
fig
```

## Figure 8. Latent heat of vaporization of water

```math
L = 2.5012 \times 10^{6} - 2378.7 \, t \qquad (-20 < t < 60)
```

[`enthalpy_of_vaporisation`](@ref) uses the regression on tabular data for temperatures above 0 °C and a regression for
sublimation from ice below.

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
htovpr = enthalpy_of_vaporisation.(tairs) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Latent heat of vaporization of water (J kg⁻¹)")
lines!(ax, ustrip.(tairs), ustrip.(u"J/kg", htovpr); linewidth = 2)
fig
```

## Figure 11. Standard atmospheric pressure

```math
p = 101325 \left[1 - 2.2569 \times 10^{-5} Z\right]^{5.2553} \qquad (-1000 < Z < 20000)
```

[`atmospheric_pressure`](@ref) uses the barometric formula with a constant lapse rate, which gives values that
differ slightly from this expression.

```@example dryair
alts = (-500:250:3000) .* u"m" # sequence of altitudes
patmos = atmospheric_pressure.(alts) # compute values

fig, ax = figure_axis("Altitude (m)", "Standard atmospheric pressure (Pa)")
lines!(ax, ustrip.(alts), ustrip.(u"Pa", patmos); linewidth = 2, label = "atmospheric_pressure")
lines!(ax, ustrip.(alts), 101325 .* (1 .- 2.2569e-5 .* ustrip.(alts)) .^ 5.2553; linewidth = 2, linestyle = :dash, label = "Smithsonian tables")
axislegend(ax; position = :rt)
fig
```

## Figure 12. Temperature

```math
t_F = \left(\frac{9}{5}\right) t_C + 32, \qquad t_C = \left(\frac{5}{9}\right)(t_F - 32)
```

Unitful converts between temperature scales:

```@example dryair
tc = collect(-20:5:50) .* u"°C" # sequence of air temperatures in °C
tf = uconvert.(u"°F", tc)
tk = uconvert.(u"K", tc)

fig, ax = figure_axis("Temperature (°C)", "Temperature (°F)")
lines!(ax, ustrip.(tc), ustrip.(tf); linewidth = 2)
fig
```

## Figure 13. Temperature coefficient of volume expansion for dry air

```math
\beta = \frac{1}{t + 273.15}
```

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
tcoeff = 1 ./ uconvert.(u"K", tairs) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Temperature coefficient of volume expansion (K⁻¹)")
lines!(ax, ustrip.(tairs), ustrip.(u"K^-1", tcoeff); linewidth = 2)
fig
```

## Figure 14. Thermal conductivity of air

```math
k = 0.02425 + 7.038 \times 10^{-5} t \qquad (-20 < t < 40)
```

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
thcond = getproperty.(dry_air_properties.(tairs, 101325.0u"Pa"), :thermal_conductivity) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Thermal conductivity of air (W m⁻¹ K⁻¹)")
lines!(ax, ustrip.(tairs), ustrip.(u"W/m/K", thcond); linewidth = 2)
fig
```

## Figure 19. Wavelength of maximum emittance from black-body

```math
\lambda_m = \frac{2.897 \times 10^{3}}{t + 273.15}
```

The wavelength is in micrometres in this expression, and [`dry_air_properties`](@ref) returns metres.

```@example dryair
tairs = collect(-20:5:50) .* u"°C" # sequence of air temperatures
emtmax = getproperty.(dry_air_properties.(tairs, 101325.0u"Pa"), :peak_wavelength) # compute values

fig, ax = figure_axis("Air temperature (°C)", "Wavelength of maximum emittance from black-body (m)")
lines!(ax, ustrip.(tairs), ustrip.(u"m", emtmax); linewidth = 2)
fig
```
