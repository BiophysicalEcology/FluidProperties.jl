const PERCENT = 100

const STULL_COEFFICIENTS = (;
    scale=0.151977,
    offset=8.313659,
    shift=1.676331,
    humidity_term=0.00391838,
    humidity_exponent=1.5,
    humidity_scale=0.023101,
    intercept=-4.686035,
)

"""
    Stull <: WetBulbMethod

Stull (2011) empirical fit to air temperature and relative humidity, for 5-99 % humidity
and -20 to 50 °C at sea level. Pressure is ignored. Use with [`wet_bulb_temperature`](@ref).

Stull, R. (2011). Wet-bulb temperature from relative humidity and air temperature.
Journal of Applied Meteorology and Climatology 50: 2267-2269.
"""
struct Stull <: WetBulbMethod end

@inline function wet_bulb_temperature(::Stull, air_temperature::Quantity, relative_humidity::Real, atmospheric_pressure::Quantity)
    (; scale, offset, shift, humidity_term, humidity_exponent, humidity_scale, intercept) = STULL_COEFFICIENTS
    # Empirical fit in °C and percent
    T = ustrip(u"°C", air_temperature)
    RH = PERCENT * _check_fraction("relative_humidity", relative_humidity)
    T_w = T * atan(scale * sqrt(RH + offset)) + atan(T + RH) - atan(RH - shift) +
        humidity_term * RH^humidity_exponent * atan(humidity_scale * RH) + intercept
    return u"K"(T_w * u"°C")
end
