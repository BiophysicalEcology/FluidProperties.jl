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

@inline function _check_fraction(name, value)
    (value < 0 || value > 1) && _throw_not_fraction(name, value)
    return value
end

# Illinois false-position root of `residual` on [T_lo, T_hi]. NaN without a sign change.
@inline function _find_root(residual, T_lo, T_hi, tolerance, max_iterations)
    f_lo = residual(T_lo)
    f_hi = residual(T_hi)
    f_hi == zero(f_hi) && return T_hi
    (f_lo > zero(f_lo)) == (f_hi > zero(f_hi)) && return T_hi * NaN
    T = T_hi
    side = 0
    for _ in 1:max_iterations
        T_previous = T
        T = (f_lo * T_hi - f_hi * T_lo) / (f_lo - f_hi)
        f = residual(T)
        (f == zero(f) || abs(T - T_previous) < tolerance) && return T
        if (f > zero(f)) == (f_lo > zero(f_lo))
            T_lo, f_lo = T, f
            side == -1 && (f_hi /= 2)
            side = -1
        else
            T_hi, f_hi = T, f
            side == 1 && (f_lo /= 2)
            side = 1
        end
    end
    return T
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
no wet bulb temperature exists between 100 K below and 1 K above the air temperature.

# Example

```julia
wet_bulb_temperature(30.0u"°C", 0.5, 101325.0u"Pa")
wet_bulb_temperature(Stull(), 30.0u"°C", 0.5, 101325.0u"Pa")
wet_bulb_temperature(Smithsonian(; vapour_pressure_equation=Huang()), 30.0u"°C", 0.5, 101325.0u"Pa")
```
"""
wet_bulb_temperature(::WetBulbMethod, ::Union{Missing,Quantity}, ::Union{Missing,Real}, ::Union{Missing,Quantity}) = missing
wet_bulb_temperature(::Union{Missing,Quantity}, ::Union{Missing,Real,SpecificHumidity}, ::Union{Missing,Quantity}; kw...) = missing
@inline wet_bulb_temperature(air_temperature::Quantity, humidity::Union{Real,SpecificHumidity}, atmospheric_pressure::Quantity; kw...) =
    wet_bulb_temperature(DaviesJones(; kw...), air_temperature, humidity, atmospheric_pressure)
