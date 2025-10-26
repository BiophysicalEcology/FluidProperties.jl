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
```julia
using Unitful

P = atmospheric_pressure(1500u"m")
"""
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
    wet_air_properties(T_drybulb, T_wetbulb, rh, T_dew, P_atmos, fO2, fCO2, fN2)
    wet_air_properties(T_drybulb; kw...)

Calculates several properties of humid air as output variables below. The program
is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
Institution Press. Washington, DC. wet_air_properties must be used in conjunction with function vapour_pressure.

Input variables are shown below. The user must supply known values for T_drybulb and P (P at one standard
atmosphere is 101 325 pascals). Values for the remaining variables are determined by whether the user has
either (1) psychrometric data (T_wetbulb or rh), or (2) hygrometric data (T_dew)

# TODO fix this desctiption
(1) Psychrometric data:
If T_wetbulb is known but not rh, then set rh=-1 and dp=999
If rh is known but not T_wetbulb then set T_wetbulb=nothing and dp=999

(2) Hygrometric data:
If T_dew is known then set T_wetublb = nothing and rh = nothing.

# Arguments

- `T_drybulb`: Dry bulb temperature (K or °C)
- `T_wetbulb`: Wet bulb temperature (K or °C)
- `rh`: Relative humidity (fractional)
- `T_dew`: Dew point temperature (K or °C)
- `P`: Barometric pressure (Pa)
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
@inline function wet_air_properties(T_drybulb; 
    T_wetbulb=nothing, 
    T_dew=nothing, 
    rh=nothing, 
    P_atmos=101325u"Pa",
    fO2=0.2095,
    fCO2=0.0004,
    fN2=0.79,
    vapour_pressure_equation=GoffGratch(),
)
    return wet_air_properties(T_drybulb, T_wetbulb, rh, T_dew, P_atmos, fO2, fCO2, fN2; vapour_pressure_equation)
end
@inline function wet_air_properties(T_drybulb, T_wetbulb, rh, T_dew, P_atmos, fO2, fCO2, fN2; vapour_pressure_equation)
    c_p_H2O_vap = 1864.40u"J/K/kg"
    c_p_dry_air = 1004.84u"J/K/kg" # should be 1006?
    f_w = 1.0053 # (-) correction factor for the departure of the mixture of air and water vapour from ideal gas laws
    M_w = (1molH₂O |> u"kg") / 1u"mol" # molar mass of water
    M_a = (fO2*molO₂ + fCO2*molCO₂ + fN2*molN₂) / 1u"mol" # molar mass of air
    P_vap_sat = vapour_pressure(vapour_pressure_equation, u"K"(T_drybulb))
    if isnothing(T_dew)
        if isnothing(rh)
            if isnothing(T_wetbulb) # We assume T_wetbulb == T_drybulb
                P_vap = P_vap_sat
                rh = 1.0
            else
                δ_bulb = u"K"(T_drybulb - T_wetbulb)
                δ_P_vap = (0.000660 * (1 + 0.00115 * ustrip(u"°C", T_wetbulb)) * ustrip(u"Pa", P) * ustrip(u"°C", δ_bulb))u"Pa"
                P_vap = vapour_pressure(vapour_pressure_equation, u"K"(T_wetbulb)) - δ_P_vap
                rh = (P_vap / P_vap_sat)
            end
        else
            P_vap = P_vap_sat * rh
        end
    else
        P_vap = vapour_pressure(vapour_pressure_equation, u"K"(T_dew))
        # And why dont we check isnothing(rh) here as well?
        rh = (P_vap / P_vap_sat)
    end
    r_w = ((M_w / M_a) * f_w * P_vap) / (P_atmos - f_w * P_vap)
    ρ_vap = P_vap * M_w / (0.998 * Unitful.R * u"K"(T_drybulb)) # TODO 0.998 a correction factor?
    ρ_vap = uconvert(u"kg/m^3", ρ_vap) # simplify units
    T_vir = u"K"(T_drybulb) * ((1 + r_w / (M_w / M_a)) / (1 + r_w))
    T_vinc = T_vir - u"K"(T_drybulb)
    ρ_air = (M_a / Unitful.R) * P_atmos / (0.999 * T_vir) # TODO 0.999 a correction factor?
    ρ_air = uconvert(u"kg/m^3", ρ_air) # simplify units
    c_p = (c_p_dry_air + (r_w * c_p_H2O_vap)) / (1 + r_w)
    ψ = if min(rh) <= 0
        -999.0u"Pa"
    else
        (4.615e+5 * ustrip(u"K", T_drybulb) * log(rh))u"Pa"
    end

    return (; P_vap, P_vap_sat, ρ_vap, r_w, T_vinc, ρ_air, c_p, ψ, rh)
end

"""
    dry_air_properties(T_drybulb; kw...)
    dry_air_properties(T_drybulb, P_atmos, elevation, fO2, fCO2, fN2)

"""
@inline dry_air_properties(T_drybulb; P_atmos=nothing, elevation=0.0u"m", fO2=0.2095, fCO2=0.0004, fN2=0.79) = 
    dry_air_properties(T_drybulb, P_atmos, elevation, fO2, fCO2, fN2)
@inline function dry_air_properties(T_drybulb, P_atmos, elevation, fO2, fCO2, fN2)
    M_a = ((fO2*molO₂ + fCO2*molCO₂ + fN2*molN₂) |> u"kg") / 1u"mol" # molar mass of air
    if isnothing(P_atmos)
        P_atmos = atmospheric_pressure(elevation)
    end
    ρ_air = (M_a / R) * P_atmos / (u"K"(T_drybulb)) # density of air
    ρ_air = uconvert(u"kg/m^3", ρ_air) # simplify units
    μ_0 = 1.8325e-5u"kg/m/s" # reference dynamic viscosity
    T_0 = 296.16u"K" # reference temperature
    C = 120u"K" # Sutherland's constant
    μ = (μ_0 * (T_0 + C) / (u"K"(T_drybulb) + C)) * (u"K"(T_drybulb) / T_0)^1.5 # dynamic viscosity, kg / m.s
    ν = μ / ρ_air # kinematic viscosity m2 / s or J.s/kg
    D_0 = 2.26e-5u"m^2/s" # reference molecular diffusivity of water vapour at 273.15 K
    D_w = D_0 * ((u"K"(T_drybulb) / 273.15u"K")^1.81) * (1.e5u"Pa" / P_atmos) # vapour diffusivity m2 / s
    k_air = (0.02425 + (7.038e-5 * ustrip(u"°C", T_drybulb)))u"W/m/K" # thermal conductivity of air
    β = 1 / u"K"(T_drybulb) # thermal expansion coefficient
    Grashof_group = g_n * β / (ν^2) # multipy by ΔT L^3 to get Grashof number, 1 / m3.K
    blackbody_emission = σ * (u"K"(T_drybulb)^4) # W/m2
    λ_max = 2.897e-3u"K*m" / (u"K"(T_drybulb)) # wavelength of maximum emission, m

    return (; P_atmos, ρ_air, μ, ν, D_w, k_air, Grashof_group, blackbody_emission, λ_max)
end

"""
    enthalpy_of_vaporisation(T::Quantity)
"""
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

    return (;
        c_p_H2O = c_p,
        ρ_H2O = ρ,
        k_H2O = K,
        μ_H2O = μ,
    )
end
