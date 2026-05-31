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
function vapour_pressure(::Huang, T)
    Tc = ustrip(u"°C", T)
    if Tc > 0.0
        # Huang (2018), water over liquid surface
        return exp((34.494 - 4924.99 / (Tc + 237.1))) / ((Tc + 105.0)^1.57) * 1u"Pa"
    else
        # Huang (2018), water over ice surface
        return exp((43.494 - 6545.8 / (Tc + 278.0))) / ((Tc + 868.0)^2) * 1u"Pa"
    end
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
function vapour_pressure(vpl::VapourPressureLookup, T)
    x  = (T - vpl.tmin) / vpl.step         # fractional 0-based index (dimensionless)
    i  = clamp(floor(Int, x) + 1, 1, length(vpl.lookup) - 1)
    w  = x - floor(x)                      # interpolation weight [0, 1)
    return vpl.lookup[i] * (1 - w) + vpl.lookup[i+1] * w
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
vapour_pressure(T) = vapour_pressure(GoffGratch(), T)
