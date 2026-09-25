const PSYCHROMETRIC_BRACKET = (; below=100.0u"K", above=1.0u"K")

# Root of `residual` between T - below and T + above. The upper margin brackets saturated air.
@inline _find_wet_bulb(residual, T, method) = _find_root(
    residual, T - PSYCHROMETRIC_BRACKET.below, T + PSYCHROMETRIC_BRACKET.above,
    method.tolerance, method.max_iterations,
)

# Temperatures in K relative to freezing, converted from the °F form (681 + 0.24 td - 0.6 tw)
const BARENBRUG_CONSTANTS = (;
    denominator=371.93u"K",
    dry_bulb_numerator=0.24,
    wet_bulb_numerator=0.6,
    dry_bulb_denominator=0.04,
    wet_bulb_denominator=0.4,
)

"""
    Barenbrug <: WetBulbMethod

    Barenbrug(; vapour_pressure_equation=GoffGratch(), tolerance=1e-5u"K", max_iterations=100)

Psychrometric equation of Barenbrug (1947, eqn 39),
`e = [eₛ(Tw) (c + 0.24 T - 0.6 Tw) - 0.24 (T - Tw) P] / (c + 0.04 T - 0.4 Tw)`
with temperatures relative to freezing, solved for the wet bulb temperature.
Use with [`wet_bulb_temperature`](@ref).

Barenbrug, A. W. T. (1947). Psychrometry and psychrometric charts. Journal of the Chemical,
Metallurgical and Mining Society of South Africa, May: 393-417.

## Keywords

- `vapour_pressure_equation`: [`VapourPressureEquation`](@ref) for saturation vapour pressure
- `tolerance`: Convergence tolerance on the temperature (K)
- `max_iterations`: Maximum number of iterations
"""
@kwdef struct Barenbrug{E,TOL} <: WetBulbMethod
    vapour_pressure_equation::E = GoffGratch()
    tolerance::TOL = 1e-5u"K"
    max_iterations::Int = 100
end

@inline function wet_bulb_temperature(method::Barenbrug, air_temperature::Quantity, relative_humidity::Real, atmospheric_pressure::Quantity)
    (; denominator, dry_bulb_numerator, wet_bulb_numerator,
       dry_bulb_denominator, wet_bulb_denominator) = BARENBRUG_CONSTANTS
    RH = _check_fraction("relative_humidity", relative_humidity)
    T = float(u"K"(air_temperature))
    P = float(u"Pa"(atmospheric_pressure))
    equation = method.vapour_pressure_equation
    T_a = T - freezing_temperature
    e = RH * vapour_pressure(equation, T)
    residual(T_w) = let T_wc = T_w - freezing_temperature
        e - (vapour_pressure(equation, T_w) * (denominator + dry_bulb_numerator * T_a - wet_bulb_numerator * T_wc) -
             dry_bulb_numerator * (T_a - T_wc) * P) /
            (denominator + dry_bulb_denominator * T_a - wet_bulb_denominator * T_wc)
    end
    return _find_wet_bulb(residual, T, method)
end

const SMITHSONIAN_CONSTANTS = (; coefficient=0.00066 / u"K", temperature_coefficient=0.00115 / u"K")

"""
    Smithsonian <: WetBulbMethod

    Smithsonian(; vapour_pressure_equation=GoffGratch(), tolerance=1e-5u"K", max_iterations=100)

Psychrometric equation of Ferrel as given in List (1971, Table 98, eqn 3),
`e = eₛ(Tw) - 0.00066 (1 + 0.00115 Tw) P (T - Tw)` with `Tw` in °C, solved for the wet bulb
temperature. Use with [`wet_bulb_temperature`](@ref).
With `vapour_pressure_equation=Bolton()` this is eqn 4 of Lemke and Kjellstrom (2012).

List, R. J. (1971). Smithsonian Meteorological Tables, 6th ed. Smithsonian Institution Press.

Lemke, B. and Kjellstrom, T. (2012). Calculating workplace WBGT from meteorological data:
a tool for climate change assessment. Industrial Health 50: 267-278.

## Keywords

- `vapour_pressure_equation`: [`VapourPressureEquation`](@ref) for saturation vapour pressure
- `tolerance`: Convergence tolerance on the temperature (K)
- `max_iterations`: Maximum number of iterations
"""
@kwdef struct Smithsonian{E,TOL} <: WetBulbMethod
    vapour_pressure_equation::E = GoffGratch()
    tolerance::TOL = 1e-5u"K"
    max_iterations::Int = 100
end

@inline function wet_bulb_temperature(method::Smithsonian, air_temperature::Quantity, relative_humidity::Real, atmospheric_pressure::Quantity)
    (; coefficient, temperature_coefficient) = SMITHSONIAN_CONSTANTS
    RH = _check_fraction("relative_humidity", relative_humidity)
    T = float(u"K"(air_temperature))
    P = float(u"Pa"(atmospheric_pressure))
    equation = method.vapour_pressure_equation
    e = RH * vapour_pressure(equation, T)
    residual(T_w) = e - (vapour_pressure(equation, T_w) -
        coefficient * (1 + temperature_coefficient * (T_w - freezing_temperature)) * P * (T - T_w))
    return _find_wet_bulb(residual, T, method)
end

"""
    EnergyBalance <: WetBulbMethod

    EnergyBalance(; vapour_pressure_equation=GoffGratch(), tolerance=1e-5u"K", max_iterations=100)

Enthalpy balance `cₚ T + L q = cₚ Tw + L qₛ(Tw)` of Zhang et al. (2021), solved for the wet bulb
temperature. `cₚ` is the specific heat of humid air, `L` the [`enthalpy_of_vaporisation`](@ref)
at the air temperature, and `q`, `qₛ` the actual and saturation specific humidities
`ε e / (P - (1 - ε) e)` with `ε` the molar mass ratio of water to air.
Use with [`wet_bulb_temperature`](@ref).

Zhang, Y., Held, I. and Fueglistaler, S. (2021). Projections of tropical heat stress constrained
by atmospheric dynamics. Nature Geoscience 14: 133-137.

## Keywords

- `vapour_pressure_equation`: [`VapourPressureEquation`](@ref) for saturation vapour pressure
- `tolerance`: Convergence tolerance on the temperature (K)
- `max_iterations`: Maximum number of iterations
"""
@kwdef struct EnergyBalance{E,TOL} <: WetBulbMethod
    vapour_pressure_equation::E = GoffGratch()
    tolerance::TOL = 1e-5u"K"
    max_iterations::Int = 100
end

@inline function wet_bulb_temperature(method::EnergyBalance, air_temperature::Quantity, relative_humidity::Real, atmospheric_pressure::Quantity)
    ε = DAVIES_JONES_CONSTANTS.molar_mass_ratio
    RH = _check_fraction("relative_humidity", relative_humidity)
    T = float(u"K"(air_temperature))
    P = float(u"Pa"(atmospheric_pressure))
    equation = method.vapour_pressure_equation
    air = wet_air_properties(T, RH, P; vapour_pressure_equation=equation)
    c_p = air.specific_heat
    specific_humidity(e) = ε * e / (P - (1 - ε) * e)
    q = specific_humidity(air.vapour_pressure)
    L = enthalpy_of_vaporisation(T)
    residual(T_w) = c_p * (T - T_w) - L * (specific_humidity(vapour_pressure(equation, T_w)) - q)
    return _find_wet_bulb(residual, T, method)
end
