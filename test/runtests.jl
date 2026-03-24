using Aqua, DataFrames, CSV, FluidProperties, Test, SafeTestsets, Unitful, UnitfulMoles

@testset "Aqua.jl quality assurance" begin
    Aqua.test_all(FluidProperties)
end

@testset "atmospheric_pressure" begin
    elevation = 100.0u"m"
    @test atmospheric_pressure(elevation) == 100128.83927387102u"Pa"
end

@safetestset "Test against NicheMapR outputs" begin include("nichemapr_tests.jl") end
@safetestset "Huang (2018) reference values" begin include("huang_tests.jl") end
