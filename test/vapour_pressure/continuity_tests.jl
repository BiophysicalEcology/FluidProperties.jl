using FluidProperties, Test, Unitful

# Ice and water branches meet where each equation switches between them
@testset "$(nameof(typeof(equation))) is continuous at freezing" for equation in (GoffGratch(), Huang())
    below = vapour_pressure(equation, 273.15u"K" - 1e-9u"K")
    above = vapour_pressure(equation, 273.15u"K" + 1e-9u"K")
    @test isapprox(below, above; rtol = 1e-8)
end
