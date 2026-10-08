const PERCENT = 100

const STULL_COEFFICIENTS = (;
    scale=0.151977,           # c₁
    offset=8.313659,          # c₂
    shift=1.676331,           # c₃
    humidity_term=0.00391838, # c₄
    humidity_exponent=1.5,    # n
    humidity_scale=0.023101,  # c₅
    intercept=-4.686035,      # c₆
)

"""
    Stull <: WetBulbMethod

Stull (2011, eqn 1) empirical fit to air temperature `T` (°C) and relative humidity `RH` (%),

`Tw = T atan(c₁ √(RH + c₂)) + atan(T + RH) - atan(RH - c₃) + c₄ RH^(3/2) atan(c₅ RH) + c₆`,

for 5-99 % humidity and -20 to 50 °C at sea level. Pressure is ignored.
Use with [`wet_bulb_temperature`](@ref).

Stull, R. (2011). Wet-bulb temperature from relative humidity and air temperature.
Journal of Applied Meteorology and Climatology 50: 2267-2269.
"""
struct Stull <: WetBulbMethod end

function wet_bulb_temperature(::Stull, air_temperature::Quantity, relative_humidity::Real, ::Quantity)
    (; scale, offset, shift, humidity_term, humidity_exponent, humidity_scale, intercept) = STULL_COEFFICIENTS

    # Empirical fit in °C and percent
    T = ustrip(u"°C", air_temperature)
    RH = PERCENT * _check_fraction("relative_humidity", relative_humidity)

    c₁ = scale
    c₂ = offset
    c₃ = shift
    c₄ = humidity_term
    c₅ = humidity_scale
    c₆ = intercept
    n = humidity_exponent

    T_w = T * atan(c₁ * sqrt(RH + c₂)) + atan(T + RH) - atan(RH - c₃) + c₄ * RH^n * atan(c₅ * RH) + c₆

    return u"K"(T_w * u"°C")
end
