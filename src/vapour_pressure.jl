abstract type VapourPressureEquation end

"""
    Teten <: VapourPressureEquation

Tetens equations for [`vapour_pressure`](@ref).

Low accuracy but very fast, with only a single `exp` call.
"""
struct Teten <: VapourPressureEquation end

vapour_pressure(::Teten, ::Missing) = missing
function vapour_pressure(::Teten, T)
    Tc = ustrip(u"°C", T)
    P_triple = 6.1071    # hPa, vapour pressure at triple point
    return P_triple * exp((17.269 * Tc) / (237.3 + Tc)) * 100u"Pa"
end

"""
    GoffGratch <: VapourPressureEquation

Widely used Goff-Gratch equations for [`vapour_pressure`](@ref).
"""
struct GoffGratch <: VapourPressureEquation end

vapour_pressure(::GoffGratch, ::Missing) = missing
function vapour_pressure(::GoffGratch, T)
    Tc = ustrip(u"°C", T)
    T = ustrip(u"K", T) + 0.01 # triple point of water is 273.16

    # Physical reference points
    T_triple = 273.16    # K, triple point of water
    T_boiling = 373.16   # K, normal boiling point of water
    P_triple = 6.1071    # hPa, vapour pressure at triple point
    P_boiling = 1013.246 # hPa, vapour pressure at boiling point

    if T < T_triple
        # --- Goff–Gratch saturation over ice ---
        logP_vap = -9.09718 * (T_triple / T - 1) +
                   -3.56654 * log10(T_triple / T) +
                   0.876793 * (1 - T / T_triple) +
                   log10(P_triple)
    else
        # --- Goff–Gratch saturation over liquid water ---
        logP_vap = -7.90298 * (T_boiling / T - 1) +
                   5.02808 * log10(T_boiling / T) +
                   -1.3816e-07 * (exp10(11.344 * (1 - T / T_boiling)) - 1) +
                   8.1328e-03 * (exp10(-3.49149 * (T_boiling / T - 1)) - 1) +
                   log10(P_boiling)
    end
    # Note: exp10 is faster than 10^x
    result = exp10(logP_vap) * 100u"Pa"

    return result
end

"""
    Huang <: VapourPressureEquation

Huang (2018) equations for [`vapour_pressure`](@ref).

High accuracy from -100 to 100 °C and reasonable performance.
"""
struct Huang <: VapourPressureEquation end

vapour_pressure(::Huang, ::Missing) = missing
function vapour_pressure(::Huang, Tk)
    t = ustrip(u"°C", Tk)
    if t > 0.0
        # Huang (2018), water over liquid surface
        return exp((34.494 - 4924.99 / (t + 237.1)) / ((t + 105.0)^1.57)) * 100u"Pa"
    else
        # Huang (2018), water over ice surface
        return exp((43.494 - 6545.8 / (t + 278.0)) / ((t + 868.0)^2)) * 100u"Pa"
    end
end

"""
    vapour_pressure(T)
    vapour_pressure(formulation, T)

Calculates saturation vapour pressure (Pa) for a given air temperature.

# Arguments
- `T`: air temperature in K.

The `GoffGratch` formulation is used by default.
"""
vapour_pressure(::Missing) = missing
vapour_pressure(Tk) = vapour_pressure(GoffGratch(), Tk)
