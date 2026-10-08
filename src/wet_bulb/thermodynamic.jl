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

# Constants of Davies-Jones (2008) and Bolton (1980), with their symbols. Mixing ratios are in kg/kg.
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
    k1_coefficients=(-53.737, 137.81, -38.5),       # k₁(π) in K, eqn 4.3
    k2_coefficients=(-0.384, 56.831, -4.392),       # k₂(π) in K, eqn 4.4
    warm_correction=1.21u"K",                       # eqns 4.10-4.11
    hot_correction=2.66u"K",                        # eqn 4.11
    hot_inverse_coefficient=0.58u"K",               # eqn 4.11
    warm_boundary=0.4,                              # (C/T_E)^λ, eqns 4.10-4.11
    minimum_equivalent_temperature=200.0u"K",
    maximum_equivalent_temperature=600.0u"K",
    newton_tolerance=0.01u"K",
    newton_iterations=4,
    newton_max_step=10.0u"K",
)

# NaN for negative bases
_nonnegative_power(x, y) = x < zero(x) ? oftype(x, NaN) : x^y

# Saturation terms along a pseudoadiabat at `temperature` τ and `pressure` p:
# f(τ, π) (Davies-Jones eqn 2.3) and f′ = ∂f/∂τ (eqns A.1-A.5)
function _pseudoadiabat_terms(temperature, pressure)
    (; heat_capacity_ratio, poisson_exponent, molar_mass_ratio, reference_pressure,
       equivalent_temperature_scale, equivalent_linear_coefficient, equivalent_quadratic_coefficient) = DAVIES_JONES_CONSTANTS

    τ = temperature
    p = pressure

    C = freezing_temperature
    λ = heat_capacity_ratio
    κ = poisson_exponent
    ε = molar_mass_ratio
    p₀ = reference_pressure
    k₀ = equivalent_temperature_scale
    k₁ = equivalent_linear_coefficient
    k₂ = equivalent_quadratic_coefficient

    (; saturation_vapour_pressure, saturation_vapour_pressure_gradient) = _vapour_pressure_terms(Bolton(), τ)
    eₛ = saturation_vapour_pressure
    deₛdτ = saturation_vapour_pressure_gradient

    π = _nonnegative_power(p / p₀, κ) # nondimensional pressure
    p₀πλ = p₀ * _nonnegative_power(π, λ)

    rₛ = ε * eₛ / (p₀πλ - eₛ + eps(Float64) * u"Pa") # saturation mixing ratio
    drₛdτ = ε * p / (p - eₛ)^2 * deₛdτ # eqn A.4

    G = (k₀ / τ - k₁) * (rₛ + k₂ * rₛ^2) # eqn 2.2
    dGdτ = -k₀ * (rₛ + k₂ * rₛ^2) / τ^2 +(k₀ / τ - k₁) * (1 + 2k₂ * rₛ) * drₛdτ # eqn A.3

    f = _nonnegative_power(C / τ, λ) * _nonnegative_power(1 - eₛ / p₀πλ, κ * λ) * exp(-λ * G) # eqn 2.3
    dlnfdτ = -λ * (1 / τ + κ * deₛdτ / (p - eₛ) + dGdτ) # eqn A.2
    f′ = f * dlnfdτ # eqn A.1

    return (;
        saturation_vapour_pressure = eₛ,
        saturation_mixing_ratio = zero(rₛ) <= rₛ <= one(rₛ) ? rₛ : oftype(rₛ, NaN),
        log_vapour_pressure_gradient = deₛdτ / eₛ,
        pseudoadiabat_function = f,
        pseudoadiabat_gradient = f′,
    )
end

# Relative humidity and mixing ratio, both fractional
function _humidity_fractions(relative_humidity::Real, saturation_mixing_ratio)
    RH = _check_fraction("relative_humidity", relative_humidity)
    rₛ = saturation_mixing_ratio
    return (; relative_humidity = RH, mixing_ratio = RH * rₛ)
end
function _humidity_fractions(specific_humidity::SpecificHumidity, saturation_mixing_ratio)
    r = _check_fraction("specific_humidity", specific_humidity.value)
    rₛ = saturation_mixing_ratio
    return (; relative_humidity = r / rₛ, mixing_ratio = r)
end

# Newton step of Davies-Jones eqn 2.6, for f(τₙ, π) = (C/T_E)^λ
function _newton_step(wet_bulb_estimate, scaled_equivalent_temperature, pressure)
    (; pseudoadiabat_function, pseudoadiabat_gradient) = _pseudoadiabat_terms(wet_bulb_estimate, pressure)
    f = pseudoadiabat_function
    f′ = pseudoadiabat_gradient
    x = scaled_equivalent_temperature
    return (f - x) / f′
end

function _solve_wet_bulb(::FixedNewton, initial_estimate, scaled_equivalent_temperature, pressure, air_temperature)
    (; newton_tolerance, newton_max_step, newton_iterations) = DAVIES_JONES_CONSTANTS

    τ = initial_estimate
    x = scaled_equivalent_temperature
    p = pressure
    Δτ_max = newton_max_step

    # Not a `while` loop from Δτ = Inf: Enzyme reverse mode doubles its derivative
    for _ in 1:newton_iterations
        Δτ = clamp(_newton_step(τ, x, p), -Δτ_max, Δτ_max)
        τ -= Δτ
        abs(Δτ) > newton_tolerance || break
    end
    return τ
end
function _solve_wet_bulb(convergence::RefinedNewton, initial_estimate, scaled_equivalent_temperature, pressure, air_temperature)
    (; tolerance, relaxation, max_iterations) = convergence

    x = scaled_equivalent_temperature
    p = pressure
    T = air_temperature
    ω = relaxation

    τ = _solve_wet_bulb(FixedNewton(), initial_estimate, x, p, T)
    for _ in 1:max_iterations
        Δτ = _newton_step(τ, x, p)
        abs(Δτ) > tolerance || return τ
        τ -= ω * Δτ
    end
    return abs(_newton_step(τ, x, p)) > tolerance ? T : τ
end

# Initial estimate of the wet bulb temperature, Davies-Jones eqns 4.3-4.11
function _initial_wet_bulb_estimate(equivalent_temperature, scaled_equivalent_temperature, pressure)
    (; reference_pressure, poisson_exponent, cold_regime_scale, regime_boundary_coefficients,
       k1_coefficients, k2_coefficients, warm_correction, hot_correction,
       hot_inverse_coefficient, warm_boundary) = DAVIES_JONES_CONSTANTS

    T_E = equivalent_temperature
    x = scaled_equivalent_temperature
    p = pressure

    C = freezing_temperature
    p₀ = reference_pressure
    κ = poisson_exponent
    A = cold_regime_scale
    c₁ = warm_correction
    c₂ = hot_correction
    c₃ = hot_inverse_coefficient

    π = _nonnegative_power(p / p₀, κ)
    D = 1 / evalpoly(p / p₀, regime_boundary_coefficients) # eqn 4.7
    k₁ = evalpoly(π, k1_coefficients) * u"K" # eqn 4.3
    k₂ = evalpoly(π, k2_coefficients) * u"K" # eqn 4.4

    return if x > D # eqn 4.8
        (; saturation_mixing_ratio, log_vapour_pressure_gradient) = _pseudoadiabat_terms(T_E, p)
        rₛ = saturation_mixing_ratio
        dlneₛdT = log_vapour_pressure_gradient
        T_E - A * rₛ / (1 + A * rₛ * dlneₛdT)
    elseif x >= 1 # eqn 4.9
        C + k₁ - k₂ * x
    elseif x >= warm_boundary # eqn 4.10
        C + (k₁ - c₁) - (k₂ - c₁) * x
    else # eqn 4.11
        C + (k₁ - c₂) - (k₂ - c₁) * x + c₃ / x
    end
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
function wet_bulb_properties(
    air_temperature::Quantity, humidity::Union{Real,SpecificHumidity}, atmospheric_pressure::Quantity;
    convergence::WetBulbConvergence=RefinedNewton(),
)
    (; heat_capacity_ratio, poisson_exponent, reference_pressure, equivalent_temperature_scale,
       equivalent_linear_coefficient, equivalent_quadratic_coefficient, moist_potential_exponent,
       lifting_temperature_offset, lifting_temperature_scale,
       minimum_equivalent_temperature, maximum_equivalent_temperature) = DAVIES_JONES_CONSTANTS

    T = float(u"K"(air_temperature))
    p = float(u"Pa"(atmospheric_pressure))

    C = freezing_temperature
    λ = heat_capacity_ratio
    κ = poisson_exponent
    p₀ = reference_pressure
    k₀ = equivalent_temperature_scale
    k₁ = equivalent_linear_coefficient
    k₂ = equivalent_quadratic_coefficient
    k₃ = moist_potential_exponent
    T₀ = lifting_temperature_offset
    T₁ = lifting_temperature_scale

    (; saturation_vapour_pressure, saturation_mixing_ratio) = _pseudoadiabat_terms(T, p)
    eₛ = saturation_vapour_pressure
    rₛ = saturation_mixing_ratio

    (; relative_humidity, mixing_ratio) = _humidity_fractions(humidity, rₛ)
    RH = relative_humidity
    r = mixing_ratio
    e = RH * eₛ

    π = _nonnegative_power(p / p₀, κ)
    T_L = 1 / (1 / (T - T₀) - log(RH) / T₁) + T₀ # Bolton eqn 22
    θ_DL = T * _nonnegative_power(p₀ / (p - e), κ) * _nonnegative_power(T / T_L, k₃ * r) # Bolton eqn 24
    θ_E = θ_DL * exp((k₀ / T_L - k₁) * r * (1 + k₂ * r)) # Bolton eqn 39
    T_E = θ_E * π
    x = _nonnegative_power(C / T_E, λ) # left side of Davies-Jones eqn 2.3

    if !(minimum_equivalent_temperature <= T_E <= maximum_equivalent_temperature)
        not_a_number = T_E * NaN
        return WetBulbProperties(;
            wet_bulb_temperature=not_a_number,
            equivalent_temperature=not_a_number,
            equivalent_potential_temperature=θ_E,
        )
    end

    τ₀ = _initial_wet_bulb_estimate(T_E, x, p)
    T_W = _solve_wet_bulb(convergence, τ₀, x, p, T)

    return WetBulbProperties(;
        wet_bulb_temperature=T_W,
        equivalent_temperature=T_E,
        equivalent_potential_temperature=θ_E,
    )
end

_davies_jones_temperature(method, air_temperature, humidity, atmospheric_pressure) =
    wet_bulb_properties(air_temperature, humidity, atmospheric_pressure; convergence=method.convergence).wet_bulb_temperature
wet_bulb_temperature(method::DaviesJones, air_temperature::Quantity, humidity::Real, atmospheric_pressure::Quantity) =
    _davies_jones_temperature(method, air_temperature, humidity, atmospheric_pressure)
wet_bulb_temperature(method::DaviesJones, air_temperature::Quantity, humidity::SpecificHumidity, atmospheric_pressure::Quantity) =
    _davies_jones_temperature(method, air_temperature, humidity, atmospheric_pressure)
wet_bulb_temperature(::DaviesJones, ::Union{Missing,Quantity}, ::SpecificHumidity, ::Union{Missing,Quantity}) = missing
