"""
    WetBulbMethod

Abstract supertype for wet bulb formulations, used by [`wet_bulb_temperature`](@ref).

Formulations: [`DaviesJones`](@ref), [`Stull`](@ref), [`Barenbrug`](@ref),
[`Smithsonian`](@ref) and [`EnergyBalance`](@ref).
"""
abstract type WetBulbMethod end

"""
    SpecificHumidity(value)

Specific humidity (kg/kg, fractional) as an input to [`wet_bulb_temperature`](@ref)
in place of relative humidity. Only [`DaviesJones`](@ref) accepts it, and treats it as
a mixing ratio.
"""
struct SpecificHumidity{V<:Real}
    value::V
end

@noinline _throw_not_fraction(name, value) = throw(DomainError(value,
    "$name must be a fraction between 0 and 1, got $value. Divide percentages by 100."))

function _check_fraction(name, value)
    (value < 0 || value > 1) && _throw_not_fraction(name, value)
    return value
end

# Root of `f` on [T_lo, T_hi], or NaN if there is no sign change or no convergence.
# Not differentiated, see `_find_temperature`.
function _bracketed_root(f, T_lo, T_hi, solver, xatol, max_iterations)
    f_lo = f(T_lo)
    f_hi = f(T_hi)
    sign(f_lo) * sign(f_hi) <= 0 || return T_hi * NaN # also catches NaN residuals
    return solve(ZeroProblem(f, (T_lo, T_hi)), solver; xatol, maxiters=max_iterations)
end

const NEWTON_REFINEMENT = (; steps=2, slope_step=1e-3) # slope_step in K

# Temperature between `lower` and `upper` where `residual` is zero.
# Roots.jl finds the root, then Newton steps refine it. Autodiff of the Newton steps gives the
# implicit derivative -(∂r/∂θ)/(∂r/∂T); two steps are needed for second derivatives.
# Units are stripped because Roots.jl allocates with Unitful quantities.
function _find_temperature(residual, lower, upper, solver, tolerance, max_iterations)
    (; steps, slope_step) = NEWTON_REFINEMENT

    f(T) = ustrip(residual(T * u"K"))

    T_lo = ustrip(u"K", lower)
    T_hi = ustrip(u"K", upper)
    xatol = ustrip(u"K", tolerance)
    h = slope_step

    T = _bracketed_root(f, T_lo, T_hi, solver, xatol, max_iterations)
    for _ in 1:steps
        ∂f∂T = (f(T + h) - f(T - h)) / 2h
        T -= f(T) / ∂f∂T
    end

    return T * u"K"
end

"""
    wet_bulb_temperature(air_temperature, humidity, atmospheric_pressure; kw...)
    wet_bulb_temperature(method, air_temperature, humidity, atmospheric_pressure)

Wet bulb temperature (K) using a [`WetBulbMethod`](@ref), by default
[`DaviesJones`](@ref). Keywords in the first form are passed to `DaviesJones`.

# Arguments

- `method`: a [`WetBulbMethod`](@ref)
- `air_temperature`: Air temperature (any temperature unit)
- `humidity`: Relative humidity (fractional, 0-1), or a [`SpecificHumidity`](@ref) for `DaviesJones`
- `atmospheric_pressure`: Barometric pressure (any pressure unit)

Humidity outside 0-1 throws a `DomainError`. Root-finding methods return `NaN` if
no wet bulb temperature exists between 100 K below and 1 K above the air temperature,
or if their solver does not converge within `max_iterations`.

# Example

```julia
wet_bulb_temperature(30.0u"°C", 0.5, 101325.0u"Pa")
wet_bulb_temperature(Stull(), 30.0u"°C", 0.5, 101325.0u"Pa")
wet_bulb_temperature(Smithsonian(; vapour_pressure_equation=Huang()), 30.0u"°C", 0.5, 101325.0u"Pa")
```
"""
wet_bulb_temperature(::WetBulbMethod, ::Union{Missing,Quantity}, ::Union{Missing,Real}, ::Union{Missing,Quantity}) = missing
wet_bulb_temperature(::Union{Missing,Quantity}, ::Union{Missing,Real,SpecificHumidity}, ::Union{Missing,Quantity}; kw...) = missing
wet_bulb_temperature(air_temperature::Quantity, humidity::Union{Real,SpecificHumidity}, atmospheric_pressure::Quantity; kw...) =
    wet_bulb_temperature(DaviesJones(; kw...), air_temperature, humidity, atmospheric_pressure)
