module FluidProperties

using Unitful, UnitfulMoles

using Unitful: ustrip, uconvert

# Define here instead of using PhysicalConstants.jl
# Its too slow due to BigFloat conversions and allocation
using Unitful: σ, R
const g_n = 9.80665u"m*s^-2"
const atm = uconvert(u"Pa", 1Unitful.atm)
const freezing_temperature = 273.15u"K"

export atmospheric_pressure, dry_air_properties, enthalpy_of_vaporisation, molar_enthalpy_of_vaporisation
export vapour_pressure, water_properties, wet_air_properties
export GasFractions, DryAirProperties, WetAirProperties, WaterProperties
export VapourPressureEquation, GoffGratch, Teten, Huang, Bolton, VapourPressureLookup
export wet_bulb_properties, wet_bulb_temperature, natural_wet_bulb_temperature
export WetBulbMethod, DaviesJones, Stull, Barenbrug, Smithsonian, EnergyBalance
export WetBulbProperties, SpecificHumidity, WetBulbConvergence, FixedNewton, RefinedNewton

@compound H2O
@compound O2
@compound CO2
@compound N2

include("vapour_pressure.jl")
include("fluid_properties.jl")
include("wet_bulb/abstract.jl")
include("wet_bulb/thermodynamic.jl")
include("wet_bulb/empirical.jl")
include("wet_bulb/psychrometric.jl")
include("wet_bulb/natural.jl")

function __init__()
    Unitful.register(FluidProperties)
end

end
