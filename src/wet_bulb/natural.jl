const NATURAL_WET_BULB_CONSTANTS = (;
    globe_excess=4.0u"K",           # Tg - T above which radiant heat is included
    globe_coefficient=0.25,
    # C(v) = 0.85 below still_wind, 1 above ventilated_wind, else intercept + slope log10(v)
    still_wind=0.03u"m/s",
    still_coefficient=0.85,
    ventilated_wind=3.0u"m/s",
    ventilated_coefficient=1.0,
    coefficient_intercept=0.96,
    coefficient_slope=0.069,
    # Radiant adjustment: 1.1 K below radiant_still_wind, -0.1 K above radiant_ventilated_wind,
    # else scale / v^exponent - offset
    radiant_still_wind=0.1u"m/s",
    radiant_still_adjustment=1.1u"K",
    radiant_ventilated_wind=1.0u"m/s",
    radiant_ventilated_adjustment=-0.1u"K",
    radiant_scale=0.1u"K",
    radiant_exponent=1.1,
    radiant_offset=0.2u"K",
)

@inline function _wind_coefficient(v)
    (; still_wind, still_coefficient, ventilated_wind, ventilated_coefficient,
       coefficient_intercept, coefficient_slope) = NATURAL_WET_BULB_CONSTANTS
    return if v < still_wind
        still_coefficient
    elseif v > ventilated_wind
        ventilated_coefficient
    else
        coefficient_intercept + coefficient_slope * log10(v / 1.0u"m/s")
    end
end

@inline function _radiant_adjustment(v)
    (; radiant_still_wind, radiant_still_adjustment, radiant_ventilated_wind,
       radiant_ventilated_adjustment, radiant_scale, radiant_exponent, radiant_offset) = NATURAL_WET_BULB_CONSTANTS
    return if v < radiant_still_wind
        radiant_still_adjustment
    elseif v > radiant_ventilated_wind
        radiant_ventilated_adjustment
    else
        radiant_scale / (v / 1.0u"m/s")^radiant_exponent - radiant_offset
    end
end

"""
    natural_wet_bulb_temperature(air_temperature, wet_bulb_temperature, wind_speed; globe_temperature=air_temperature)

Natural wet bulb temperature (K) from the psychrometric wet bulb temperature (for example from
[`wet_bulb_temperature`](@ref)) and wind speed `v`, after Bernard and Pourmoghani (1999, eqns 1 and 2).

Without radiant heat (globe temperature `Tg` less than 4 K above air temperature `T`):

`Tnwb = T - C (T - Tw)`, with `C = 0.85` for `v < 0.03` m/s, `C = 1` for `v > 3` m/s,
and `C = 0.96 + 0.069 log10(v)` between.

With radiant heat:

`Tnwb = Tw + 0.25 (Tg - T) + e`, with `e = 1.1` K for `v < 0.1` m/s, `e = -0.1` K for `v > 1` m/s,
and `e = 0.1 / v^1.1 - 0.2` K between.

The default `globe_temperature` is the air temperature. Wind speed must not be negative.

Bernard, T. E. and Pourmoghani, M. (1999). Prediction of workplace wet bulb global temperature.
Applied Occupational and Environmental Hygiene 14: 126-134.
"""
natural_wet_bulb_temperature(::Union{Missing,Quantity}, ::Union{Missing,Quantity}, ::Union{Missing,Quantity}; kw...) = missing
@inline function natural_wet_bulb_temperature(
    air_temperature::Quantity, wet_bulb_temperature::Quantity, wind_speed::Quantity;
    globe_temperature::Quantity=air_temperature,
)
    T = float(u"K"(air_temperature))
    T_w = float(u"K"(wet_bulb_temperature))
    T_g = float(u"K"(globe_temperature))
    v = float(u"m/s"(wind_speed))
    v < zero(v) && throw(DomainError(v, "wind_speed must not be negative, got $v"))
    return if T_g - T >= NATURAL_WET_BULB_CONSTANTS.globe_excess
        T_w + NATURAL_WET_BULB_CONSTANTS.globe_coefficient * (T_g - T) + _radiant_adjustment(v)
    else
        T - _wind_coefficient(v) * (T - T_w)
    end
end
