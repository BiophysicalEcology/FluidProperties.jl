module FluidProperties

using Unitful, UnitfulMoles

using Unitful: ustrip, uconvert

# Define here instead of using PhysicalConstants.jl
# Its too slow due to BigFloat conversions and allocation
using Unitful: σ, R
const g_n = 9.80665u"m*s^-2"
const atm = uconvert(u"Pa", 1Unitful.atm)

export atmospheric_pressure, dry_air_properties, enthalpy_of_vaporisation, molar_enthalpy_of_vaporisation
export vapour_pressure, water_properties, wet_air_properties
export GasFractions, DryAirProperties, WetAirProperties, WaterProperties

@compound H2O
@compound O2
@compound CO2
@compound N2

include("vapour_pressure.jl")
include("fluid_properties.jl")

function __init__()\
    Unitful.register(FluidProperties)
end

end
