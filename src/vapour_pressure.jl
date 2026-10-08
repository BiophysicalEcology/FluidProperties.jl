"""
    VapourPressureEquation

Abstract supertype for the saturation vapour pressure equations used by [`vapour_pressure`](@ref):
[`GoffGratch`](@ref), [`Teten`](@ref), [`Huang`](@ref), [`Bolton`](@ref), [`ClausiusClapeyron`](@ref)
and [`VapourPressureLookup`](@ref).
"""
abstract type VapourPressureEquation end

Base.broadcastable(equation::VapourPressureEquation) = Ref(equation)

"""
    Teten <: VapourPressureEquation

Tetens equations for [`vapour_pressure`](@ref).

Low accuracy but very fast, with only a single `exp` call.
"""
struct Teten <: VapourPressureEquation end

vapour_pressure(::Teten, ::Missing) = missing
function vapour_pressure(::Teten, temperature)
    T_c = ustrip(u"°C", temperature)

    e₀ = 6.1071 # hPa, vapour pressure at triple point
    a = 17.269
    b = 237.3   # °C

    eₛ = e₀ * exp(a * T_c / (b + T_c)) * 100u"Pa"

    return eₛ
end

"""
    GoffGratch <: VapourPressureEquation

Widely used Goff-Gratch equations for [`vapour_pressure`](@ref).
"""
struct GoffGratch <: VapourPressureEquation end

vapour_pressure(::GoffGratch, ::Missing) = missing
function vapour_pressure(::GoffGratch, temperature)
    # Clamped to avoid solvers taking this near/below zero and causing negative logs below.
    T = clamp(ustrip(u"K", temperature) + 0.01, 173.15, 373.16) # triple point of water is 273.16

    # Physical reference points
    T_tr = 273.16    # K, triple point of water
    T_b = 373.16     # K, normal boiling point of water
    e_tr = 6.1071    # hPa, vapour pressure at triple point
    e_b = 1013.246   # hPa, vapour pressure at boiling point

    log₁₀eₛ = if T < T_tr
        # Goff–Gratch saturation over ice
        -9.09718 * (T_tr / T - 1) +
        -3.56654 * log10(T_tr / T) +
        0.876793 * (1 - T / T_tr) +
        log10(e_tr)
    else
        # Goff–Gratch saturation over liquid water
        -7.90298 * (T_b / T - 1) +
        5.02808 * log10(T_b / T) +
        -1.3816e-07 * (exp10(11.344 * (1 - T / T_b)) - 1) +
        8.1328e-03 * (exp10(-3.49149 * (T_b / T - 1)) - 1) +
        log10(e_b)
    end
    # Note: exp10 is faster than 10^x
    eₛ = exp10(log₁₀eₛ) * 100u"Pa"

    return eₛ
end

"""
    Huang <: VapourPressureEquation

Huang (2018) equations for [`vapour_pressure`](@ref).

High accuracy from -100 to 100 °C and reasonable performance.
"""
struct Huang <: VapourPressureEquation end

vapour_pressure(::Huang, ::Missing) = missing
function vapour_pressure(::Huang, temperature)
    t = ustrip(u"°C", temperature)

    eₛ = if t > 0.0
        # Huang (2018), water over liquid surface
        exp(34.494 - 4924.99 / (t + 237.1)) / (t + 105.0)^1.57 * 1u"Pa"
    else
        # Huang (2018), water over ice surface
        exp(43.494 - 6545.8 / (t + 278.0)) / (t + 868.0)^2 * 1u"Pa"
    end

    return eₛ
end

"""
    Bolton <: VapourPressureEquation

Bolton (1980) eqn 10 for [`vapour_pressure`](@ref), over liquid water.
Used by [`DaviesJones`](@ref).
"""
struct Bolton <: VapourPressureEquation end

# Saturation vapour pressure (Pa) and its temperature gradient (Pa/K), Bolton eqn 10
function _vapour_pressure_terms(::Bolton, temperature)
    T = u"K"(temperature)

    e₀ = 611.2u"Pa"
    a = 17.67
    b = 243.5u"K"
    C = freezing_temperature

    eₛ = e₀ * exp(a * (T - C) / (T - C + b))
    deₛdT = eₛ * a * b / (T - C + b)^2

    return (; saturation_vapour_pressure = eₛ, saturation_vapour_pressure_gradient = deₛdT)
end

vapour_pressure(::Bolton, ::Missing) = missing
vapour_pressure(model::Bolton, temperature) = _vapour_pressure_terms(model, temperature).saturation_vapour_pressure

"""
    ClausiusClapeyron <: VapourPressureEquation

    ClausiusClapeyron(; ice=true, kw...)

The Clausius-Clapeyron relation integrated from the triple point with constant specific heats
(the Rankine-Kirchhoff approximation), for [`vapour_pressure`](@ref):

```math
e^* = p_{tr} \\left(\\frac{T}{T_{tr}}\\right)^{Δc_p/R_v}
      \\exp\\left[\\frac{L_0 - Δc_p T_0}{R_v}\\left(\\frac{1}{T_{tr}} - \\frac{1}{T}\\right)\\right]
```

where `L₀` is the latent heat of vaporisation (or sublimation) at `T₀` and `Δcₚ` is the
difference between the specific heats of water vapour and liquid water (or ice).

This is the formulation of [Thermodynamics.jl](https://github.com/CliMA/Thermodynamics.jl),
and the defaults are the values of [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl),
so the results match `saturation_vapor_pressure(param_set, T, Liquid())` above the triple
point and `saturation_vapor_pressure(param_set, T, Ice())` below it.

## Keywords

- `ice`: Use saturation over ice below the triple point (default: `true`),
  otherwise over (supercooled) liquid water at all temperatures
- `triple_point_pressure`: Triple point pressure of water (Pa)
- `triple_point_temperature`: Triple point temperature of water (K)
- `reference_temperature`: Temperature `T₀` of the reference latent heats (K)
- `latent_heat_vaporisation`: Latent heat of vaporisation at `T₀` (J/kg)
- `latent_heat_sublimation`: Latent heat of sublimation at `T₀` (J/kg)
- `specific_heat_vapour`: Isobaric specific heat of water vapour (J/kg/K)
- `specific_heat_liquid`: Isobaric specific heat of liquid water (J/kg/K)
- `specific_heat_ice`: Isobaric specific heat of ice (J/kg/K)
- `gas_constant_vapour`: Specific gas constant of water vapour (J/kg/K)
"""
@kwdef struct ClausiusClapeyron{P,T,L,C,R} <: VapourPressureEquation
    ice::Bool = true
    triple_point_pressure::P = 611.657u"Pa"
    triple_point_temperature::T = 273.16u"K"
    reference_temperature::T = 273.16u"K"
    latent_heat_vaporisation::L = 2.5008e6u"J/kg"
    latent_heat_sublimation::L = 2.8344e6u"J/kg"
    specific_heat_vapour::C = 1859.0u"J/kg/K"
    specific_heat_liquid::C = 4181.0u"J/kg/K"
    specific_heat_ice::C = 2070.0u"J/kg/K"
    gas_constant_vapour::R = 461.5u"J/kg/K"
end

vapour_pressure(::ClausiusClapeyron, ::Missing) = missing
function vapour_pressure(cc::ClausiusClapeyron, T)
    T = u"K"(T)
    p_tr = cc.triple_point_pressure
    T_tr = cc.triple_point_temperature
    T_0 = cc.reference_temperature
    R_v = cc.gas_constant_vapour
    over_ice = cc.ice && T < T_tr
    L_0 = over_ice ? cc.latent_heat_sublimation : cc.latent_heat_vaporisation
    Δc_p = cc.specific_heat_vapour - (over_ice ? cc.specific_heat_ice : cc.specific_heat_liquid)
    return p_tr * exp(Δc_p / R_v * log(T / T_tr) + (L_0 - Δc_p * T_0) / R_v * (1 / T_tr - 1 / T))
end

"""
    VapourPressureLookup <: VapourPressureEquation

Lookup-table with linear interpolation for [`vapour_pressure`](@ref).

Pre-computes a table of saturation vapour pressures over a temperature range at
construction time. Evaluation uses linear interpolation with no transcendental
function calls, making it significantly faster than equation-based formulations
for repeated calls.

# Keyword arguments
- `Tmin`: minimum temperature (°C), default -40.0
- `Tmax`: maximum temperature (°C), default 60.0
- `dT`: temperature step size (°C), default 0.1
- `formulation`: [`VapourPressureEquation`](@ref) used to build the table, default `GoffGratch()`

# Example
```julia
vpl = VapourPressureLookup()
vpl = VapourPressureLookup(tmin=-80.0u"°C", tmax=100.0u"°C", step=0.05u"K", formulation=Huang())
es  = vapour_pressure(vpl, 293.15u"K")
```
"""
struct VapourPressureLookup<: VapourPressureEquation
    tmin::typeof(1.0u"K")
    step::typeof(1.0u"K")
    lookup::Vector{typeof(1.0u"Pa")}
end
function VapourPressureLookup(formulation=GoffGratch(); tmin=-40.0u"°C", tmax=90.0u"°C", step=0.1u"K")
    unit(step) == u"K" || throw(ArgumentError("step must be given in Kelvin (K), got $(unit(step))"))
    ts = tmin:step:tmax
    table = [vapour_pressure(formulation, t) for t in ts]
    return VapourPressureLookup(tmin, step, table)
end

vapour_pressure(::VapourPressureLookup, ::Missing) = missing
function vapour_pressure(vpl::VapourPressureLookup, temperature)
    T = temperature
    T_min = vpl.tmin
    ΔT = vpl.step
    eₛ = vpl.lookup

    x = (T - T_min) / ΔT # fractional 0-based index (dimensionless)
    i = clamp(floor(Int, x) + 1, 1, length(eₛ) - 1)
    w = x - floor(x)     # interpolation weight [0, 1)

    return eₛ[i] * (1 - w) + eₛ[i+1] * w
end

"""
    vapour_pressure(temperature)
    vapour_pressure(formulation, temperature)

Calculates saturation vapour pressure (Pa) for a given air temperature.

# Arguments
- `temperature`: air temperature (any temperature unit).

The `GoffGratch` formulation is used by default.
"""
vapour_pressure(::Missing) = missing
vapour_pressure(temperature) = vapour_pressure(GoffGratch(), temperature)
