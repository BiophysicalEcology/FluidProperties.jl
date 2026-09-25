"""
    WetBulbConvergence

Abstract supertype for the iteration strategies of [`DaviesJones`](@ref).
"""
abstract type WetBulbConvergence end

"""
    FixedNewton <: WetBulbConvergence

Newton-Raphson for at most four iterations, stopping when a step is below 0.01 K
(Davies-Jones 2008, eqn 2.6).
"""
struct FixedNewton <: WetBulbConvergence end

"""
    RefinedNewton <: WetBulbConvergence

    RefinedNewton(; tolerance=1e-5u"K", relaxation=1.0, max_iterations=20000)

[`FixedNewton`](@ref) followed by Newton-Raphson steps until the step is below
`tolerance`. Returns the air temperature if `max_iterations` is reached.

## Keywords

- `tolerance`: Convergence tolerance on the step (K)
- `relaxation`: Fraction of each step applied
- `max_iterations`: Maximum number of steps
"""
@kwdef struct RefinedNewton{TOL,RL} <: WetBulbConvergence
    tolerance::TOL = 1e-5u"K"
    relaxation::RL = 1.0
    max_iterations::Int = 20000
end

"""
    DaviesJones <: WetBulbMethod

    DaviesJones(; convergence=RefinedNewton())

Davies-Jones (2008) wet bulb method for [`wet_bulb_temperature`](@ref) and
[`wet_bulb_properties`](@ref). Accepts relative humidity or [`SpecificHumidity`](@ref).

## Keywords

- `convergence`: [`RefinedNewton`](@ref) (default) or [`FixedNewton`](@ref)
"""
@kwdef struct DaviesJones{C<:WetBulbConvergence} <: WetBulbMethod
    convergence::C = RefinedNewton()
end

"""
    WetBulbProperties

    WetBulbProperties(; wet_bulb_temperature, equivalent_temperature, equivalent_potential_temperature)

Returned by [`wet_bulb_properties`](@ref).

## Fields

- `wet_bulb_temperature`: Wet bulb temperature (K)
- `equivalent_temperature`: Equivalent temperature (K)
- `equivalent_potential_temperature`: Equivalent potential temperature (K)
"""
@kwdef struct WetBulbProperties{W,E,P}
    wet_bulb_temperature::W
    equivalent_temperature::E
    equivalent_potential_temperature::P
end

Base.show(io::IO, ::MIME"text/plain", properties::WetBulbProperties) = _show_properties(io, properties)
Base.show(io::IO, properties::WetBulbProperties) = _show_properties_compact(io, properties)

# Symbols are those of Davies-Jones (2008) and Bolton (1980). Mixing ratios are in kg/kg.
const DAVIES_JONES_CONSTANTS = (;
    heat_capacity_ratio=3.504,                      # λ = cₚd/R_d
    poisson_exponent=0.2854,                        # κ
    molar_mass_ratio=0.6220,                        # ε, water over dry air
    reference_pressure=100000.0u"Pa",               # p₀
    equivalent_temperature_scale=3036.0u"K",        # k₀, Bolton eqn 39
    equivalent_linear_coefficient=1.78,             # k₁
    equivalent_quadratic_coefficient=0.448,         # k₂
    moist_potential_exponent=0.28,                  # k₃, Bolton eqn 24
    lifting_temperature_offset=55.0u"K",            # Bolton eqn 22
    lifting_temperature_scale=2840.0u"K",           # Bolton eqn 22
    cold_regime_scale=2675.0u"K",                   # A, eqn 4.8
    regime_boundary_coefficients=(0.6512, 0.1859),  # D(p) = 1 / (c₀ + c₁ p/p₀), eqn 4.7
    k1_coefficients=(-53.737, 137.81, -38.5),       # k₁(Π) in K, eqn 4.3
    k2_coefficients=(-0.384, 56.831, -4.392),       # k₂(Π) in K, eqn 4.4
    warm_correction=1.21u"K",                       # eqns 4.10-4.11
    hot_correction=2.66u"K",                        # eqn 4.11
    hot_inverse_coefficient=0.58u"K",               # eqn 4.11
    warm_boundary=0.4,                              # x, eqns 4.10-4.11
    valid_equivalent_temperature=(200.0u"K", 600.0u"K"),
    newton_tolerance=0.01u"K",
    newton_iterations=4,
    newton_max_step=10.0u"K",
)

# NaN for negative bases
@inline _pow(x, y) = x < zero(x) ? oftype(x, NaN) : x^y

# Saturation terms at temperature `T` and pressure `p`: eₛ, rₛ (kg/kg, NaN outside 0-1),
# dln(eₛ)/dT, f(T) (Davies-Jones eqn 2.3) and f′ = ∂f/∂T (eqn A.1)
@inline function _saturation_terms(T, p)
    C = freezing_temperature
    λ = DAVIES_JONES_CONSTANTS.heat_capacity_ratio
    κ = DAVIES_JONES_CONSTANTS.poisson_exponent
    ε = DAVIES_JONES_CONSTANTS.molar_mass_ratio
    p₀ = DAVIES_JONES_CONSTANTS.reference_pressure
    k₀ = DAVIES_JONES_CONSTANTS.equivalent_temperature_scale
    k₁ = DAVIES_JONES_CONSTANTS.equivalent_linear_coefficient
    k₂ = DAVIES_JONES_CONSTANTS.equivalent_quadratic_coefficient

    eₛ, deₛ = _vapour_pressure_terms(Bolton(), T)
    dlneₛ = deₛ / eₛ
    p_e = p - eₛ

    Π = _pow(p / p₀, κ)
    p₀Πλ = p₀ * _pow(Π, λ)
    r = ε * eₛ / (p₀Πλ - eₛ + eps(Float64) * u"Pa")
    drdT = ε * p / (p_e * p_e) * deₛ

    G = (k₀ / T - k₁) * (r + k₂ * r * r)
    dGdT = -k₀ * (r + k₂ * r * r) / (T * T) + (k₀ / T - k₁) * (1 + 2k₂ * r) * drdT

    f = _pow(C / T, λ) * _pow(1 - eₛ / p₀Πλ, κ * λ) * exp(-λ * G)
    dlnfdT = -λ * (1 / T + κ * deₛ / p_e + dGdT)
    f′ = f * dlnfdT

    rₛ = zero(r) <= r <= one(r) ? r : oftype(r, NaN)
    return (; eₛ, rₛ, dlneₛ, f, f′)
end

# Relative humidity RH and mixing ratio r, both fractional
@inline function _humidity_fractions(relative_humidity::Real, rₛ)
    RH = _check_fraction("relative_humidity", relative_humidity)
    return RH, RH * rₛ
end
@inline function _humidity_fractions(specific_humidity::SpecificHumidity, rₛ)
    r = _check_fraction("specific_humidity", specific_humidity.value)
    return r / rₛ, r
end

# Newton step for f(T_W) = x, Davies-Jones eqn 2.6
@inline function _newton_step(T_W, x, p)
    (; f, f′) = _saturation_terms(T_W, p)
    return (f - x) / f′
end

@inline function _solve_wet_bulb(::FixedNewton, T_W, x, p, T)
    tolerance = DAVIES_JONES_CONSTANTS.newton_tolerance
    max_step = DAVIES_JONES_CONSTANTS.newton_max_step
    ΔT = oftype(T_W, Inf * u"K")
    i = 0
    while abs(ΔT) > tolerance && i < DAVIES_JONES_CONSTANTS.newton_iterations
        ΔT = clamp(_newton_step(T_W, x, p), -max_step, max_step)
        T_W -= ΔT
        i += 1
    end
    return T_W
end
@inline function _solve_wet_bulb(convergence::RefinedNewton, T_W, x, p, T)
    T_W = _solve_wet_bulb(FixedNewton(), T_W, x, p, T)
    (; tolerance, relaxation, max_iterations) = convergence
    ΔT = _newton_step(T_W, x, p)
    i = 0
    while abs(ΔT) > tolerance && i < max_iterations
        T_W -= relaxation * ΔT
        ΔT = _newton_step(T_W, x, p)
        i += 1
    end
    return abs(ΔT) > tolerance ? T : T_W
end

"""
    wet_bulb_properties(air_temperature, humidity, atmospheric_pressure; convergence=RefinedNewton())

Wet bulb, equivalent and equivalent potential temperatures after Davies-Jones (2008),
from the equivalent potential temperature of Bolton (1980, eqns 22, 24, 39).
Saturation vapour pressure is [`Bolton`](@ref).

# Arguments

- `air_temperature`: Air temperature (any temperature unit)
- `humidity`: Relative humidity (fractional, 0-1) or a [`SpecificHumidity`](@ref)
- `atmospheric_pressure`: Barometric pressure (any pressure unit)

# Keywords

- `convergence`: [`RefinedNewton`](@ref) (default) or [`FixedNewton`](@ref)

Humidity outside 0-1 throws a `DomainError`.

# Returns

A [`WetBulbProperties`](@ref). The wet bulb and equivalent temperatures are `NaN`
where the equivalent temperature is outside 200-600 K.

# Example

```julia
wet_bulb_properties(30.0u"°C", 0.5, 101325.0u"Pa")
wet_bulb_properties(30.0u"°C", SpecificHumidity(0.0135), 101325.0u"Pa"; convergence=FixedNewton())
```

# References

- Bolton, D. (1980). The computation of equivalent potential temperature.
  Monthly Weather Review 108(7): 1046-1053.
- Davies-Jones, R. (2008). An efficient and accurate method for computing the
  wet-bulb temperature along pseudoadiabats. Monthly Weather Review 136(7): 2764-2785.
"""
wet_bulb_properties(::Union{Missing,Quantity}, ::Union{Missing,Real,SpecificHumidity}, ::Union{Missing,Quantity}; kw...) = missing
@inline function wet_bulb_properties(
    air_temperature::Quantity, humidity::Union{Real,SpecificHumidity}, atmospheric_pressure::Quantity;
    convergence::WetBulbConvergence=RefinedNewton(),
)
    constants = DAVIES_JONES_CONSTANTS
    C = freezing_temperature
    λ = constants.heat_capacity_ratio
    κ = constants.poisson_exponent
    p₀ = constants.reference_pressure
    k₀ = constants.equivalent_temperature_scale
    k₁ = constants.equivalent_linear_coefficient
    k₂ = constants.equivalent_quadratic_coefficient
    k₃ = constants.moist_potential_exponent
    A = constants.cold_regime_scale
    T_offset = constants.lifting_temperature_offset
    T_scale = constants.lifting_temperature_scale
    T_E_min, T_E_max = constants.valid_equivalent_temperature

    T = float(u"K"(air_temperature))
    p = float(u"Pa"(atmospheric_pressure))

    (; eₛ, rₛ) = _saturation_terms(T, p)
    RH, r = _humidity_fractions(humidity, rₛ)
    e = eₛ * RH

    # Regime parameters, eqns 4.3, 4.4, 4.7
    Π = _pow(p / p₀, κ)
    D = 1 / evalpoly(p / p₀, constants.regime_boundary_coefficients)
    kπ₁ = evalpoly(Π, constants.k1_coefficients) * u"K"
    kπ₂ = evalpoly(Π, constants.k2_coefficients) * u"K"

    T_L = 1 / (1 / (T - T_offset) - log(RH) / T_scale) + T_offset     # Bolton eqn 22
    θ_DL = T * _pow(p₀ / (p - e), κ) * _pow(T / T_L, k₃ * r)          # Bolton eqn 24
    θ_E = θ_DL * exp((k₀ / T_L - k₁) * r * (1 + k₂ * r))              # Bolton eqn 39
    T_E = θ_E * Π
    x = _pow(C / T_E, λ)

    if !(T_E_min <= T_E <= T_E_max)
        not_a_number = T_E * NaN
        return WetBulbProperties(;
            wet_bulb_temperature=not_a_number,
            equivalent_temperature=not_a_number,
            equivalent_potential_temperature=θ_E,
        )
    end

    # First guess of T_W, eqns 4.8-4.11
    T_W = if x > D
        cold = _saturation_terms(T_E, p)
        T_E - A * cold.rₛ / (1 + A * cold.rₛ * cold.dlneₛ)
    elseif x >= 1
        C + kπ₁ - kπ₂ * x
    elseif x >= constants.warm_boundary
        C + (kπ₁ - constants.warm_correction) - (kπ₂ - constants.warm_correction) * x
    else
        C + (kπ₁ - constants.hot_correction) - (kπ₂ - constants.warm_correction) * x +
            constants.hot_inverse_coefficient / x
    end

    T_W = _solve_wet_bulb(convergence, T_W, x, p, T)

    return WetBulbProperties(;
        wet_bulb_temperature=T_W,
        equivalent_temperature=T_E,
        equivalent_potential_temperature=θ_E,
    )
end

@inline _davies_jones_temperature(method, air_temperature, humidity, atmospheric_pressure) =
    wet_bulb_properties(air_temperature, humidity, atmospheric_pressure; convergence=method.convergence).wet_bulb_temperature
@inline wet_bulb_temperature(method::DaviesJones, air_temperature::Quantity, humidity::Real, atmospheric_pressure::Quantity) =
    _davies_jones_temperature(method, air_temperature, humidity, atmospheric_pressure)
@inline wet_bulb_temperature(method::DaviesJones, air_temperature::Quantity, humidity::SpecificHumidity, atmospheric_pressure::Quantity) =
    _davies_jones_temperature(method, air_temperature, humidity, atmospheric_pressure)
wet_bulb_temperature(::DaviesJones, ::Union{Missing,Quantity}, ::SpecificHumidity, ::Union{Missing,Quantity}) = missing
