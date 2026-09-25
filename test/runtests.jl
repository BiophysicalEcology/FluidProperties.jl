using Aqua, FluidProperties, Test, SafeTestsets, Unitful, UnitfulMoles

@testset "Aqua.jl quality assurance" begin
    Aqua.test_all(FluidProperties)
end

@safetestset "Fluid properties" begin include("fluid_properties_tests.jl") end
@safetestset "Huang (2018) reference values" begin include("vapour_pressure/huang_tests.jl") end
@safetestset "VapourPressureLookup" begin include("vapour_pressure/lookup_tests.jl") end
@safetestset "Bolton vapour pressure" begin include("vapour_pressure/bolton_tests.jl") end
@safetestset "Clausius-Clapeyron vapour pressure" begin include("vapour_pressure/clausius_clapeyron_tests.jl") end
@safetestset "Davies-Jones wet bulb" begin include("wet_bulb/thermodynamic_tests.jl") end
@safetestset "Wet bulb methods" begin include("wet_bulb/methods_tests.jl") end
@safetestset "Natural wet bulb" begin include("wet_bulb/natural_tests.jl") end
