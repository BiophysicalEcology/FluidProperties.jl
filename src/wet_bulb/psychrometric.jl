const PSYCHROMETRIC_BRACKET = (; below=100.0u"K", above=1.0u"K")

# Root of `residual` between T - below and T + above. The upper margin brackets saturated air.
function _find_wet_bulb(residual, air_temperature, method)
    (; below, above) = PSYCHROMETRIC_BRACKET
    (; solver, tolerance, max_iterations) = method
    T = air_temperature
    return _find_temperature(residual, T - below, T + above, solver, tolerance, max_iterations)
end

# Temperatures in K relative to freezing, converted from the °F form (681 + 0.24 td - 0.6 tw)
const BARENBRUG_CONSTANTS = (;
    denominator=371.93u"K",     # c₀
    dry_bulb_numerator=0.24,    # c₁
    wet_bulb_numerator=0.6,     # c₂
    dry_bulb_denominator=0.04,  # c₃
    wet_bulb_denominator=0.4,   # c₄
)

"""
    Barenbrug <: WetBulbMethod

    Barenbrug(; vapour_pressure_equation=GoffGratch(), solver=A42(), tolerance=1e-5u"K", max_iterations=100)

Psychrometric equation of Barenbrug (1947, eqn 39),
`e = [eₛ(Tw) (c₀ + c₁ T - c₂ Tw) - c₁ (T - Tw) P] / (c₀ + c₃ T - c₄ Tw)`
with temperatures relative to freezing, `c₀ = 371.93` K, `c₁ = 0.24`, `c₂ = 0.6`,
`c₃ = 0.04` and `c₄ = 0.4`, solved for the wet bulb temperature.
Use with [`wet_bulb_temperature`](@ref).

Barenbrug, A. W. T. (1947). Psychrometry and psychrometric charts. Journal of the Chemical,
Metallurgical and Mining Society of South Africa, May: 393-417.

## Keywords

- `vapour_pressure_equation`: [`VapourPressureEquation`](@ref) for saturation vapour pressure
- `solver`: Roots.jl bracketing method: `A42()` (default), `AlefeldPotraShi()`, `FalsePosition()` or `Bisection()`
- `tolerance`: Convergence tolerance on the temperature (K)
- `max_iterations`: Maximum number of iterations
"""
@kwdef struct Barenbrug{E,S,TOL} <: WetBulbMethod
    vapour_pressure_equation::E = GoffGratch()
    solver::S = A42()
    tolerance::TOL = 1e-5u"K"
    max_iterations::Int = 100
end

function wet_bulb_temperature(method::Barenbrug, air_temperature::Quantity, relative_humidity::Real, atmospheric_pressure::Quantity)
    (; denominator, dry_bulb_numerator, wet_bulb_numerator,
       dry_bulb_denominator, wet_bulb_denominator) = BARENBRUG_CONSTANTS

    T = float(u"K"(air_temperature))
    RH = _check_fraction("relative_humidity", relative_humidity)
    P = float(u"Pa"(atmospheric_pressure))
    equation = method.vapour_pressure_equation

    C = freezing_temperature
    c₀ = denominator
    c₁ = dry_bulb_numerator
    c₂ = wet_bulb_numerator
    c₃ = dry_bulb_denominator
    c₄ = wet_bulb_denominator

    t = T - C
    e = RH * vapour_pressure(equation, T)
    residual(T_w) = let t_w = T_w - C, eₛ = vapour_pressure(equation, T_w)
        e - (eₛ * (c₀ + c₁ * t - c₂ * t_w) - c₁ * (t - t_w) * P) / (c₀ + c₃ * t - c₄ * t_w)
    end

    return _find_wet_bulb(residual, T, method)
end

const SMITHSONIAN_CONSTANTS = (;
    psychrometer_coefficient=0.00066 / u"K",  # A
    temperature_coefficient=0.00115 / u"K",   # B
)

"""
    Smithsonian <: WetBulbMethod

    Smithsonian(; vapour_pressure_equation=GoffGratch(), solver=A42(), tolerance=1e-5u"K", max_iterations=100)

Psychrometric equation of Ferrel as given in List (1971, Table 98, eqn 3),
`e = eₛ(Tw) - A (1 + B Tw) P (T - Tw)` with `Tw` in °C, `A = 0.00066` K⁻¹ and
`B = 0.00115` K⁻¹, solved for the wet bulb temperature. Use with [`wet_bulb_temperature`](@ref).
With `vapour_pressure_equation=Bolton()` this is eqn 4 of Lemke and Kjellstrom (2012).

List, R. J. (1971). Smithsonian Meteorological Tables, 6th ed. Smithsonian Institution Press.

Lemke, B. and Kjellstrom, T. (2012). Calculating workplace WBGT from meteorological data:
a tool for climate change assessment. Industrial Health 50: 267-278.

## Keywords

- `vapour_pressure_equation`: [`VapourPressureEquation`](@ref) for saturation vapour pressure
- `solver`: Roots.jl bracketing method: `A42()` (default), `AlefeldPotraShi()`, `FalsePosition()` or `Bisection()`
- `tolerance`: Convergence tolerance on the temperature (K)
- `max_iterations`: Maximum number of iterations
"""
@kwdef struct Smithsonian{E,S,TOL} <: WetBulbMethod
    vapour_pressure_equation::E = GoffGratch()
    solver::S = A42()
    tolerance::TOL = 1e-5u"K"
    max_iterations::Int = 100
end

function wet_bulb_temperature(method::Smithsonian, air_temperature::Quantity, relative_humidity::Real, atmospheric_pressure::Quantity)
    (; psychrometer_coefficient, temperature_coefficient) = SMITHSONIAN_CONSTANTS

    T = float(u"K"(air_temperature))
    RH = _check_fraction("relative_humidity", relative_humidity)
    P = float(u"Pa"(atmospheric_pressure))
    equation = method.vapour_pressure_equation

    C = freezing_temperature
    A = psychrometer_coefficient
    B = temperature_coefficient

    e = RH * vapour_pressure(equation, T)
    residual(T_w) = let eₛ = vapour_pressure(equation, T_w)
        e - (eₛ - A * (1 + B * (T_w - C)) * P * (T - T_w))
    end

    return _find_wet_bulb(residual, T, method)
end

"""
    EnergyBalance <: WetBulbMethod

    EnergyBalance(; vapour_pressure_equation=GoffGratch(), solver=A42(), tolerance=1e-5u"K", max_iterations=100)

Enthalpy balance `cₚ T + L q = cₚ Tw + L qₛ(Tw)` of Zhang et al. (2021), solved for the wet bulb
temperature. `cₚ` is the specific heat of humid air, `L` the [`enthalpy_of_vaporisation`](@ref)
at the air temperature, and `q`, `qₛ` the actual and saturation specific humidities
`ε e / (P - (1 - ε) e)` with `ε` the molar mass ratio of water to air.
Use with [`wet_bulb_temperature`](@ref).

Zhang, Y., Held, I. and Fueglistaler, S. (2021). Projections of tropical heat stress constrained
by atmospheric dynamics. Nature Geoscience 14: 133-137.

## Keywords

- `vapour_pressure_equation`: [`VapourPressureEquation`](@ref) for saturation vapour pressure
- `solver`: Roots.jl bracketing method: `A42()` (default), `AlefeldPotraShi()`, `FalsePosition()` or `Bisection()`
- `tolerance`: Convergence tolerance on the temperature (K)
- `max_iterations`: Maximum number of iterations
"""
@kwdef struct EnergyBalance{E,S,TOL} <: WetBulbMethod
    vapour_pressure_equation::E = GoffGratch()
    solver::S = A42()
    tolerance::TOL = 1e-5u"K"
    max_iterations::Int = 100
end

function wet_bulb_temperature(method::EnergyBalance, air_temperature::Quantity, relative_humidity::Real, atmospheric_pressure::Quantity)
    T = float(u"K"(air_temperature))
    RH = _check_fraction("relative_humidity", relative_humidity)
    P = float(u"Pa"(atmospheric_pressure))
    equation = method.vapour_pressure_equation

    ε = DAVIES_JONES_CONSTANTS.molar_mass_ratio

    air = wet_air_properties(T, RH, P; vapour_pressure_equation=equation)
    c_p = air.specific_heat
    e = air.vapour_pressure
    L = enthalpy_of_vaporisation(T)

    q = ε * e / (P - (1 - ε) * e)
    residual(T_w) = let eₛ = vapour_pressure(equation, T_w), qₛ = ε * eₛ / (P - (1 - ε) * eₛ)
        c_p * (T - T_w) - L * (qₛ - q)
    end

    return _find_wet_bulb(residual, T, method)
end
