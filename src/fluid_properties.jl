"""
    GasFractions

    GasFractions(fO2, fCO2, fN2)
    GasFractions(; fO2, fCO2, fN2)

Atmospheric gas composition as mole fractions.
Default values represent standard dry air.

## Fields / Keywords

- `fO2`: Oxygen fraction (default: 0.2095)
- `fCO2`: Carbon dioxide fraction (default: 0.0004)
- `fN2`: Nitrogen fraction (default: 0.79)
"""
Base.@kwdef struct GasFractions{O,C,N}
    fO2::O = 0.2095
    fCO2::C = 0.0004
    fN2::N = 0.79
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
    atmospheric_pressure(h::Quantity;
                 h_ref::Quantity = 0u"m",
                 P_ref::Quantity = 101325u"Pa",
                 L_ref::Quantity = -0.0065u"K/m",
                 T_ref::Quantity = 288u"K",
                 g_0::Quantity = 9.80665u"m/s^2",
                 M::Quantity = 0.0289644u"kg/mol") -> Quantity

Computes atmospheric pressure at a given altitude using the barometric formula,
assuming a constant temperature lapse rate (standard tropospheric approximation).

# Arguments
- `h`: Elevation at which to compute pressure (with length units, e.g. `u"m"`).
- `h_ref`: Reference altitude (default: `0u"m"`).
- `P_ref`: Pressure at the reference altitude (default: `101325u"Pa"`, standard sea-level pressure).
- `L_ref`: Temperature lapse rate (default: `-0.0065u"K/m"`).
- `T_ref`: Temperature at the reference altitude (default: `288u"K"`).
- `g_0`: Standard gravitational acceleration (default: `9.80665u"m/s^2"`).
- `M`: Molar mass of dry air (default: `0.0289644u"kg/mol"`).

# Returns
- `P_a`: Atmospheric pressure at altitude `h` (with pressure units, e.g. `u"Pa"`).

# Notes
- Uses the simplified barometric formula assuming a linear lapse rate and ideal gas behavior.
- The universal gas constant `R` is used from the `Unitful` package.

# Example
using Unitful

P = atmospheric_pressure(1500u"m")
"""
atmospheric_pressure(::Missing; kw...) = missing
function atmospheric_pressure(h::Quantity;
    #P_ref::Quantity = atm,
    L_ref::Quantity = -0.0065u"K/m",
    T_ref::Quantity = 288.0u"K",
    #g_0::Quantity = g_n,
    M::Quantity = 0.0289644u"kg/mol"
)
    P_a = atm * (1 + (L_ref / T_ref) * h) ^ ((-g_n * M) / (R * L_ref)) #5.2553026003237262u"kg*m^2*J^-1*s^-2"#

    return P_a
end

"""
    wet_air_properties(T, rh, P; gasfrac=GasFractions(), vapour_pressure_equation=GoffGratch())
    wet_air_properties(T_drybulb; kw...)

Calculates several properties of humid air as output variables below. The program
is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
Institution Press. Washington, DC. wet_air_properties must be used in conjunction with function vapour_pressure.

Input variables are shown below. The user must supply known values for T_drybulb and P (P at one standard
atmosphere is 101 325 pascals). Values for the remaining variables are determined by whether the user has
either (1) psychrometric data (T_wetbulb or rh), or (2) hygrometric data (T_dew)

# Arguments

- `T_drybulb`: Dry bulb temperature (K or °C)
- `rh`: Relative humidity (fractional)
- `P_atmos`: Barometric pressure (Pa)
- `fO2`; fractional O2 concentration in atmosphere, -
- `fCO2`; fractional CO2 concentration in atmosphere, -
- `fN2`; fractional N2 concentration in atmosphere, -
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
@inline wet_air_properties(T, rh, P;
    gasfrac::GasFractions=GasFractions(),
    vapour_pressure_equation=GoffGratch(),
) = wet_air_properties(T, rh, P, gasfrac, vapour_pressure_equation)
@inline function wet_air_properties(
    T_drybulb::Quantity,
    rh::Real,
    P_atmos::Quantity,
    gasfrac::GasFractions,
    vapour_pressure_equation,
)
    (; fO2, fCO2, fN2) = gasfrac
    T = u"K"(T_drybulb)
    P = P_atmos

    # Constants
    c_p_H2O_vap = 1864.40u"J/K/kg"
    c_p_dry_air = 1004.84u"J/K/kg"
    f_w = 1.0053  # correction factor for departure from ideal gas laws

    # Molar masses
    M_w = u"kg"(1molH₂O) / 1u"mol"  # water
    M_a = (fO2 * molO₂ + fCO2 * molCO₂ + fN2 * molN₂) / 1u"mol"  # air

    # Vapour pressure
    P_sat = FluidProperties.vapour_pressure(vapour_pressure_equation, T)
    P_v = P_sat * rh

    # Mixing ratio
    r = ((M_w / M_a) * f_w * P_v) / (P - f_w * P_v)
    mixing_ratio = r

    # Vapour density
    ρ_v = P_v * M_w / (0.998 * Unitful.R * T)
    vapour_density = uconvert(u"kg/m^3", ρ_v)

    # Virtual temperature and increment
    T_vir = T * ((1 + r / (M_w / M_a)) / (1 + r))
    virtual_temp_increment = T_vir - T

    # Air density
    ρ = (M_a / Unitful.R) * P / (0.999 * T_vir)
    density = uconvert(u"kg/m^3", ρ)

    # Specific heat
    specific_heat = (c_p_dry_air + (r * c_p_H2O_vap)) / (1 + r)

    # Water potential
    water_potential = rh <= 0 ? -999.0u"Pa" : (4.615e+5 * ustrip(u"K", T) * log(rh))u"Pa"


    return WetAirProperties(;
        density,
        specific_heat,
        vapour_pressure=P_v,
        vapour_density,
        mixing_ratio,
        relative_humidity=rh,
        water_potential,
        virtual_temp_increment,
    )
end

"""
    dry_air_properties(T, P; gasfrac=GasFractions())
"""
dry_air_properties(::Missing, ::Missing; kwargs...) = missing
dry_air_properties(::Missing; kwargs...) = missing

@inline dry_air_properties(T, P; gasfrac::GasFractions=GasFractions()) =
    dry_air_properties(T, P, gasfrac)

@inline dry_air_properties(T; P_atmos=101325u"Pa", gasfrac::GasFractions=GasFractions()) =
    dry_air_properties(T, P_atmos, gasfrac)

@inline function dry_air_properties(T_drybulb::Quantity, P_atmos::Quantity, gasfrac::GasFractions)
    T = u"K"(T_drybulb)
    P = P_atmos

    # Molar mass of air
    (; fO2, fCO2, fN2) = gasfrac
    M_a = (fO2 * molO₂ + fCO2 * molCO₂ + fN2 * molN₂) / 1u"mol"

    # Density
    ρ = uconvert(u"kg/m^3", (M_a / R) * P / T)

    # Dynamic viscosity (Sutherland's formula)
    μ_0 = 1.8325e-5u"kg/m/s"  # reference dynamic viscosity
    T_0 = 296.16u"K"          # reference temperature
    C = 120.0u"K"             # Sutherland's constant
    μ = (μ_0 * (T_0 + C) / (T + C)) * (T / T_0)^1.5

    # Kinematic viscosity
    ν = μ / ρ

    # Thermal conductivity
    thermal_conductivity = (0.02425 + (7.038e-5 * ustrip(u"°C", T_drybulb)))u"W/m/K"

    # Vapour diffusivity
    D_0 = 2.26e-5u"m^2/s"  # reference at 273.15 K
    vapour_diffusivity = D_0 * ((T / 273.15u"K")^1.81) * (1.e5u"Pa" / P)

    # Grashof coefficient (multiply by ΔT·L³ to get Grashof number)
    β = 1 / T  # thermal expansion coefficient
    grashof_coefficient = g_n * β / (ν^2)

    # Blackbody radiation
    blackbody_emission = σ * T^4
    peak_wavelength = 2.897e-3u"K*m" / T

    return DryAirProperties(;
        density=ρ,
        molar_mass=M_a,
        dynamic_viscosity=μ,
        kinematic_viscosity=ν,
        thermal_conductivity,
        vapour_diffusivity,
        grashof_coefficient,
        blackbody_emission,
        peak_wavelength,
    )
end

"""
    enthalpy_of_vaporisation(T::Quantity)
"""
enthalpy_of_vaporisation(::Missing) = missing
function enthalpy_of_vaporisation(T::Quantity)
    # These regressions don't respect units, so we strip them
    # convert any temperature (K or °C) to Celsius
    T = ustrip(u"°C", T)
    if T > 0
        return u"J/kg"((2500.8 - 2.36 * T + 0.0016 * T^2 - 0.00006 * T^3) * u"kJ/kg")
    else
        return u"J/kg"((834.1 - 0.29 * T - 0.004 * T^2) * u"kJ/kg")
    end
end


"""
    molar_enthalpy_of_vaporisation(T::Quantity)

From Campbell et al 1994 p. 309

References
- Campbell, G. S., Jungbauer, J. D. Jr., Bidlake, W. R., & Hungerford, R. D. (1994). 
  Predicting the effect of temperature on soil thermal conductivity. 
  Soil Science, 158(5), 307–313.
"""
molar_enthalpy_of_vaporisation(::Missing) = missing
function molar_enthalpy_of_vaporisation(T::Quantity)
    # This regressions doesn't respect units, so we strip them
    # convert any temperature (K or °C) to Celsius
    T = ustrip(u"°C", T)
    return (45144.0 - 48.0 * T) * u"J/mol"
end


"""
    water_properties(T::Quantity)

Compute phyiscal properties of liquid water at a given temperature `T`.

# Description
These properties are based on regressions obtained using Grapher from Golden Software 
on data from:

> Ede, "An Introduction of Heat Transfer Principles and Calculations," Pergamon Press, 1967, p. 262.  
> (Regression performed by W. Porter, 14 July 1988)

The regressions are valid for temperatures up to 60°C. Temperatures above 60°C are clamped to 60°C.

# Inputs
- `T::Quantity`: Temperature of water. Can be specified with units (e.g., `u"°C"`). Internally converted to °C.

# Returns
A named tuple with the following fields (all returned as Unitful quantities):

- `c_p_H2O` : Specific heat capacity of water, J/(kg·K)
- `ρ_H2O`   : Density of water, kg/m^3
- `k_H2O`   : Thermal conductivity of water, W/(m·K)
- `μ_H2O`   : Dynamic viscosity of water, kg/(m·s)

# Notes
- Input temperatures are converted to °C and then units stripped before calculation.
- Density is clamped at 60°C to avoid extrapolation beyond the regression range.
"""
water_properties(::Missing) = missing
function water_properties(T::Quantity)
    # These regressions don't respect units, so we strip them
    T = ustrip(u"°C", T) # Ensure temperature is in °C

    # Specific heat capacity (J/kg·K)
    c_p = (4220.02 - 4.5531 * T + 0.182958 * T^2 - 0.00310614 * T^3 + 1.89399e-5 * T^4) * u"J/kg/K"

    # Density (kg/m^3)
    ρ = if T < 30
        1000.0 * u"kg/m^3"
    elseif T <= 60
        (1017.0 - 0.6 * T) * u"kg/m^3"
    else
        T = 60.0
        (1017.0 - 0.6 * T) * u"kg/m^3" # Clamp to 60°C
    end

    # Thermal conductivity (W/m·K)
    K = (0.551666 + 0.00282144 * T - 2.02383e-5 * T^2) * u"W/m/K"

    # Dynamic viscosity (kg/m·s)
    μ = (0.0017515 - 4.31502e-5 * T + 3.71431e-7 * T^2) * u"kg/m/s"

    return WaterProperties(;
        density=ρ,
        specific_heat=c_p,
        thermal_conductivity=K,
        dynamic_viscosity=μ,
    )
end
