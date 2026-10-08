"""
    GasFractions

    GasFractions(oxygen, carbon_dioxide, nitrogen)
    GasFractions(; oxygen, carbon_dioxide, nitrogen)

Atmospheric gas composition as mole fractions.
Default values represent standard dry air.

## Fields / Keywords

- `oxygen`: Oxygen fraction (default: 0.2095)
- `carbon`: Carbon dioxide fraction (default: 0.0004)
- `nitrogen`: Nitrogen fraction (default: 0.79)
"""
Base.@kwdef struct GasFractions{O,C,N}
    oxygen::O = 0.2095
    carbon_dioxide::C = 0.0004
    nitrogen::N = 0.79
end

"""
    DryAirProperties

    DryAirProperties(; density, molar_mass, dynamic_viscosity, ...)

Properties of dry air at a given temperature and pressure.
Returned by [`dry_air_properties`](@ref).

## Fields

- `density`: Air density (kg/m³)
- `molar_mass`: Molar mass of air (kg/mol)
- `dynamic_viscosity`: Dynamic viscosity (kg/m/s)
- `kinematic_viscosity`: Kinematic viscosity (m²/s)
- `thermal_conductivity`: Thermal conductivity (W/m/K)
- `vapour_diffusivity`: Water vapour diffusivity in air (m²/s)
- `grashof_coefficient`: Grashof coefficient, multiply by ΔT·L³ for Grashof number (1/m³/K)
- `blackbody_emission`: Blackbody radiation at temperature (W/m²)
- `peak_wavelength`: Wien's law peak emission wavelength (m)
"""
@kwdef struct DryAirProperties{D,M,DV,KV,TC,VD,GR,BB,PW}
    density::D
    molar_mass::M
    dynamic_viscosity::DV
    kinematic_viscosity::KV
    thermal_conductivity::TC
    vapour_diffusivity::VD
    grashof_coefficient::GR
    blackbody_emission::BB
    peak_wavelength::PW
end

"""
    WaterProperties

    WaterProperties(; density, specific_heat, thermal_conductivity, dynamic_viscosity)

Physical properties of liquid water at a given temperature.
Returned by [`water_properties`](@ref).

## Fields

- `density`: Water density (kg/m³)
- `specific_heat`: Specific heat capacity (J/kg/K)
- `thermal_conductivity`: Thermal conductivity (W/m/K)
- `dynamic_viscosity`: Dynamic viscosity (kg/m/s)
"""
@kwdef struct WaterProperties{D,SH,TC,DV}
    density::D
    specific_heat::SH
    thermal_conductivity::TC
    dynamic_viscosity::DV
end

"""
    WetAirProperties

    WetAirProperties(; density, specific_heat, vapour_pressure, ...)

Properties of humid air at a given temperature, relative humidity, and pressure.
Returned by [`wet_air_properties`](@ref).

## Fields

- `density`: Air density (kg/m³)
- `specific_heat`: Specific heat at constant pressure (J/kg/K)
- `vapour_pressure`: Partial pressure of water vapour (Pa)
- `vapour_density`: Water vapour density (kg/m³)
- `mixing_ratio`: Mass of water vapour per mass of dry air (kg/kg)
- `relative_humidity`: Relative humidity (fractional, 0-1)
- `water_potential`: Water potential (Pa)
- `virtual_temp_increment`: Virtual temperature minus actual temperature (K)
"""
@kwdef struct WetAirProperties{D,SH,VP,VD,MR,RH,WP,VTI}
    density::D
    specific_heat::SH
    vapour_pressure::VP
    vapour_density::VD
    mixing_ratio::MR
    relative_humidity::RH
    water_potential::WP
    virtual_temp_increment::VTI
end

"""
    atmospheric_pressure(elevation::Quantity;
                 reference_pressure::Quantity = atm,
                 laps_rate::Quantity = -0.0065u"K/m",
                 temperature::Quantity = 288u"K",
                 M::Quantity = 0.0289644u"kg/mol") -> Quantity

Computes atmospheric pressure at a given altitude using the barometric formula,
assuming a constant temperature lapse rate (standard tropospheric approximation).

# Arguments
- `elevation`: Elevation at which to compute pressure (with length units, e.g. `u"m"`).
- `reference_pressure`: Pressure at `elevation = 0` (default: standard atmosphere,
  `atm`).
- `laps_rate`: Temperature lapse rate (default: `-0.0065u"K/m"`).
- `reference_temperature`: Temperature at the altitude (default: `288u"K"`).
- `air_molar_mass`: Molar mass of dry air (default: `0.0289644u"kg/mol"`).

# Returns
Atmospheric pressure at altitude `h` (with pressure units, e.g. `u"Pa"`).

# Notes
- Uses the simplified barometric formula assuming a linear lapse rate and ideal gas behavior.
- The universal gas constant `R` is used from the `Unitful` package.
"""
atmospheric_pressure(::Missing; kw...) = missing
function atmospheric_pressure(elevation::Quantity;
    reference_pressure::Quantity = atm,
    laps_rate::Quantity = -0.0065u"K/m",
    temperature::Quantity = 288.0u"K",
    air_molar_mass::Quantity = 0.0289644u"kg/mol"
)
    z = elevation
    P₀ = reference_pressure
    Γ = laps_rate
    T₀ = temperature
    M = air_molar_mass

    P = P₀ * (1 + (Γ / T₀) * z)^((-g_n * M) / (R * Γ))

    return P
end

# Constants of the humid air equations of List (1971), as in NicheMapR WETAIR (Tracy et al. 2016)
const WET_AIR_CONSTANTS = (;
    specific_heat_dry_air=1004.84u"J/K/kg",         # c_p_a
    specific_heat_water_vapour=1864.40u"J/K/kg",    # c_p_v
    enhancement_factor=1.0053,                      # f_w, departure of humid air from the ideal gas laws
    vapour_compressibility=0.998,                   # Z_v, compressibility factor of water vapour
    air_compressibility=0.999,                      # Z_a, compressibility factor of humid air
    gas_constant_water_vapour=461.5u"J/K/kg",       # R_v
    liquid_water_density=1000.0u"kg/m^3",           # ρ_w, for water potential
    dry_air_water_potential=-999.0u"Pa",            # ψ returned where relative humidity is zero
)

"""
    wet_air_properties(T, rh, P; gas_fractions=GasFractions(), vapour_pressure_equation=GoffGratch())
    wet_air_properties(T_drybulb; kw...)

Calculates several properties of humid air as output variables below. The program
is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
Institution Press. Washington, DC. wet_air_properties must be used in conjunction with function vapour_pressure.

Input variables are shown below. The user must supply known values for T_drybulb and P (P at one standard
atmosphere is 101 325 pascals). Values for the remaining variables are determined by whether the user has
either (1) psychrometric data (T_wetbulb or rh), or (2) hygrometric data (T_dew)

# Arguments

- `drybulb_temperature`: Dry bulb temperature (K or °C)
- `relative_humidity`: Relative humidity (fractional)
- `atmospheric_pressure`: Barometric pressure (Pa)
- `oxygen`; fractional O2 concentration in atmosphere, -
- `carbon_dioxide`; fractional CO2 concentration in atmosphere, -
- `nitrogen`; fractional N2 concentration in atmosphere, -

# - `P_vap`: Vapour pressure (Pa)
# - `P_vap_sat`: Saturation vapour pressure (Pa)
# - `ρ_vap`: Vapour density (kg m-3)
# - `r_w Mixing`: ratio (kg kg-1)
# - `T_vir`: Virtual temperature (K)
# - `T_vinc`: Virtual temperature increment (K)
# - `ρ_air`: Density of the air (kg m-3)
# - `c_p`: Specific heat of air at constant pressure (J kg-1 K-1)
# - `ψ`: Water potential (Pa)
# - `rh`: Relative humidity (fractional)

"""
wet_air_properties(::Missing, ::Missing, ::Missing; kwargs...) = missing

wet_air_properties(::Missing; kwargs...) = missing
function wet_air_properties(T, rh, P;
    gas_fractions::GasFractions=GasFractions(),
    vapour_pressure_equation=GoffGratch(),
)
    wet_air_properties(T, rh, P, gas_fractions, vapour_pressure_equation)
end
function wet_air_properties(
    drybulb_temperature::Quantity,
    relative_humidity::Real,
    atmospheric_pressure::Quantity,
    gas_fractions::GasFractions,
    vapour_pressure_equation,
)
    (; specific_heat_dry_air, specific_heat_water_vapour, enhancement_factor, vapour_compressibility,
       air_compressibility, gas_constant_water_vapour, liquid_water_density, dry_air_water_potential) = WET_AIR_CONSTANTS

    T = u"K"(drybulb_temperature)
    rh = relative_humidity
    p = atmospheric_pressure

    f_O₂ = gas_fractions.oxygen
    f_CO₂ = gas_fractions.carbon_dioxide
    f_N₂ = gas_fractions.nitrogen

    c_p_a = specific_heat_dry_air
    c_p_v = specific_heat_water_vapour
    f_w = enhancement_factor
    Z_v = vapour_compressibility
    Z_a = air_compressibility
    R_v = gas_constant_water_vapour
    ρ_w = liquid_water_density
    ψ_dry = dry_air_water_potential

    # Molecular weights
    M_v = u"kg"(1molH₂O) / 1u"mol"                                # water
    M_a = (f_O₂ * molO₂ + f_CO₂ * molCO₂ + f_N₂ * molN₂) / 1u"mol" # dry air

    # Vapour pressure
    e_s = vapour_pressure(vapour_pressure_equation, T)
    e = e_s * rh

    # Mixing ratio
    r_w = ((M_v / M_a) * f_w * e) / (p - f_w * e)

    # Vapour density
    ρ_v = uconvert(u"kg/m^3", e * M_v / (Z_v * R * T))

    # Virtual temperature and increment
    T_v = T * ((1 + r_w / (M_v / M_a)) / (1 + r_w))
    ΔT_v = T_v - T

    # Air density
    ρ = uconvert(u"kg/m^3", (M_a / R) * p / (Z_a * T_v))

    # Specific heat
    c_p = (c_p_a + (r_w * c_p_v)) / (1 + r_w)

    # Water potential
    ψ = rh <= 0 ? ψ_dry : uconvert(u"Pa", ρ_w * R_v * T * log(rh))

    return WetAirProperties(;
        density=ρ,
        specific_heat=c_p,
        vapour_pressure=e,
        vapour_density=ρ_v,
        mixing_ratio=r_w,
        relative_humidity=rh,
        water_potential=ψ,
        virtual_temp_increment=ΔT_v,
    )
end

# Constants of the dry air equations, as in NicheMapR DRYAIR (Tracy et al. 2016)
const DRY_AIR_CONSTANTS = (;
    reference_viscosity=1.8325e-5u"kg/m/s",                     # μ₀, Sutherland's formula
    viscosity_reference_temperature=296.16u"K",                 # T₀
    sutherland_constant=120.0u"K",                              # C
    sutherland_exponent=1.5,                                    # m
    thermal_conductivity_coefficients=(0.02425u"W/m/K", 7.038e-5u"W/m/K^2"), # k = a + b t, t relative to freezing
    reference_diffusivity=2.26e-5u"m^2/s",                      # D₀
    diffusivity_reference_temperature=273.15u"K",               # T_D₀
    diffusivity_reference_pressure=1.0e5u"Pa",                  # p₀
    diffusivity_exponent=1.81,                                  # n
    wien_constant=2.897e-3u"K*m",                               # b, Wien's displacement law
)

"""
    dry_air_properties(T, P; gas_fractions=GasFractions())
"""
dry_air_properties(::Missing, ::Missing; kwargs...) = missing
dry_air_properties(::Missing; kwargs...) = missing
dry_air_properties(T, P; gas_fractions::GasFractions=GasFractions()) =
    dry_air_properties(T, P, gas_fractions)
dry_air_properties(T; atmospheric_pressure=atm, gas_fractions::GasFractions=GasFractions()) =
    dry_air_properties(T, atmospheric_pressure, gas_fractions)
function dry_air_properties(
    drybulb_temperature::Quantity, atmospheric_pressure::Quantity, gas_fractions::GasFractions
)
    (; reference_viscosity, viscosity_reference_temperature, sutherland_constant, sutherland_exponent,
       thermal_conductivity_coefficients, reference_diffusivity, diffusivity_reference_temperature,
       diffusivity_reference_pressure, diffusivity_exponent, wien_constant) = DRY_AIR_CONSTANTS

    T = u"K"(drybulb_temperature)
    t = T - freezing_temperature
    p = atmospheric_pressure

    f_O₂ = gas_fractions.oxygen
    f_CO₂ = gas_fractions.carbon_dioxide
    f_N₂ = gas_fractions.nitrogen

    μ₀ = reference_viscosity
    T₀ = viscosity_reference_temperature
    C = sutherland_constant
    m = sutherland_exponent
    D₀ = reference_diffusivity
    T_D₀ = diffusivity_reference_temperature
    p₀ = diffusivity_reference_pressure
    n = diffusivity_exponent
    b = wien_constant

    # Molecular weight of dry air
    M_a = (f_O₂ * molO₂ + f_CO₂ * molCO₂ + f_N₂ * molN₂) / 1u"mol"

    # Density
    ρ = uconvert(u"kg/m^3", (M_a / R) * p / T)

    # Dynamic viscosity (Sutherland's formula)
    μ = (μ₀ * (T₀ + C) / (T + C)) * (T / T₀)^m

    # Kinematic viscosity
    ν = μ / ρ

    # Thermal conductivity
    k = evalpoly(t, thermal_conductivity_coefficients)

    # Diffusivity of water vapour in air
    D = D₀ * ((T / T_D₀)^n) * (p₀ / p)

    # Group of variables in the Grashof number (multiply by ΔT·L³ to get the Grashof number)
    β = 1 / T # temperature coefficient of volume expansion
    γ = g_n * β / (ν^2)

    # Black-body emittance
    φ = σ * T^4

    # Wavelength of maximum emittance (Wien's displacement law)
    λ_m = b / T

    return DryAirProperties(;
        density=ρ,
        molar_mass=M_a,
        dynamic_viscosity=μ,
        kinematic_viscosity=ν,
        thermal_conductivity=k,
        vapour_diffusivity=D,
        grashof_coefficient=γ,
        blackbody_emission=φ,
        peak_wavelength=λ_m,
    )
end

# Polynomial coefficients in temperature relative to freezing, lowest order first
const ENTHALPY_OF_VAPORISATION_COEFFICIENTS = (;
    water=(2500.8u"kJ/kg", -2.36u"kJ/kg/K", 0.0016u"kJ/kg/K^2", -0.00006u"kJ/kg/K^3"), # above freezing
    ice=(2834.1u"kJ/kg", -0.29u"kJ/kg/K", -0.004u"kJ/kg/K^2"),                         # at and below freezing
)

"""
    enthalpy_of_vaporisation(temperature::Quantity)
"""
enthalpy_of_vaporisation(::Missing) = missing
function enthalpy_of_vaporisation(temperature::Quantity)
    (; water, ice) = ENTHALPY_OF_VAPORISATION_COEFFICIENTS

    T = u"K"(temperature)
    t = T - freezing_temperature

    L = if T > freezing_temperature
        evalpoly(t, water)
    else
        evalpoly(t, ice)
    end

    return u"J/kg"(L)
end

# Polynomial coefficients in temperature relative to freezing, lowest order first
const MOLAR_ENTHALPY_OF_VAPORISATION_COEFFICIENTS = (45144.0u"J/mol", -48.0u"J/mol/K")

"""
    molar_enthalpy_of_vaporisation(T::Quantity)

From Campbell et al 1994 p. 309

References
- Campbell, G. S., Jungbauer, J. D. Jr., Bidlake, W. R., & Hungerford, R. D. (1994).
  Predicting the effect of temperature on soil thermal conductivity.
  Soil Science, 158(5), 307–313.
"""
molar_enthalpy_of_vaporisation(::Missing) = missing
function molar_enthalpy_of_vaporisation(temperature::Quantity)
    T = u"K"(temperature)
    t = T - freezing_temperature

    λ = evalpoly(t, MOLAR_ENTHALPY_OF_VAPORISATION_COEFFICIENTS)

    return λ
end

# Regressions of Porter (1988) on Ede (1967). Polynomial coefficients in temperature
# relative to freezing, lowest order first.
const WATER_PROPERTIES_CONSTANTS = (;
    specific_heat_coefficients=(4220.02u"J/kg/K", -4.5531u"J/kg/K^2", 0.182958u"J/kg/K^3",
                                -0.00310614u"J/kg/K^4", 1.89399e-5u"J/kg/K^5"),
    cold_density=1000.0u"kg/m^3",                              # below density_threshold
    density_coefficients=(1017.0u"kg/m^3", -0.6u"kg/m^3/K"),   # from density_threshold to maximum_temperature
    density_threshold=303.15u"K",                              # 30 °C
    maximum_temperature=333.15u"K",                            # 60 °C, upper limit of the regressions
    thermal_conductivity_coefficients=(0.551666u"W/m/K", 0.00282144u"W/m/K^2", -2.02383e-5u"W/m/K^3"),
    dynamic_viscosity_coefficients=(0.0017515u"kg/m/s", -4.31502e-5u"kg/m/s/K", 3.71431e-7u"kg/m/s/K^2"),
)

"""
    water_properties(T::Quantity)

Compute phyiscal properties of liquid water at a given temperature `T`.

# Description
These properties are based on regressions obtained using Grapher from Golden Software
on data from:

> Ede, "An Introduction of Heat Transfer Principles and Calculations," Pergamon Press, 1967, p. 262.
> (Regression performed by W. Porter, 14 July 1988)

The regressions are valid for temperatures up to 60°C. Above 60°C, density, thermal
conductivity and dynamic viscosity are those at 60°C.

# Inputs
- `T::Quantity`: Temperature of water. Can be specified with units (e.g., `u"°C"`).

# Returns

A `WaterProperties` object, with the following fields (all Unitful quantities):

- `density`   : Density of water, kg/m^3
- `specific_heat` : Specific heat capacity of water, J/(kg·K)
- `thermal_conductivity`   : Thermal conductivity of water, W/(m·K)
- `dynamic_viscosity`   : Dynamic viscosity of water, kg/(m·s)
"""
water_properties(::Missing) = missing
function water_properties(temperature::Quantity)
    (; specific_heat_coefficients, cold_density, density_coefficients, density_threshold,
       maximum_temperature, thermal_conductivity_coefficients, dynamic_viscosity_coefficients) = WATER_PROPERTIES_CONSTANTS

    T = u"K"(temperature)

    T_ρ = density_threshold
    T_max = maximum_temperature
    ρ_cold = cold_density

    t = T - freezing_temperature
    t_clamped = min(T, T_max) - freezing_temperature # clamped to the upper limit of the regressions

    # Specific heat capacity
    c_p = evalpoly(t, specific_heat_coefficients)

    # Density
    ρ = T < T_ρ ? ρ_cold : evalpoly(t_clamped, density_coefficients)

    # Thermal conductivity
    k = evalpoly(t_clamped, thermal_conductivity_coefficients)

    # Dynamic viscosity
    μ = evalpoly(t_clamped, dynamic_viscosity_coefficients)

    return WaterProperties(;
        density=ρ,
        specific_heat=c_p,
        thermal_conductivity=k,
        dynamic_viscosity=μ,
    )
end
