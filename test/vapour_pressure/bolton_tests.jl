using FluidProperties
using Test
using Unitful

@testset "Bolton vapour pressure" begin
    @test vapour_pressure(Bolton(), 20.0u"°C") ≈ 2338.0u"Pa" rtol=1e-3
    @test unit(vapour_pressure(Bolton(), 293.15u"K")) == u"Pa"
    @test ismissing(vapour_pressure(Bolton(), missing))
    # Gradient against a finite difference
    _, deₛ = FluidProperties._vapour_pressure_terms(Bolton(), 20.0u"°C")
    h = 1e-5u"K"
    T = 293.15u"K"
    fd = (vapour_pressure(Bolton(), T + h) - vapour_pressure(Bolton(), T - h)) / 2h
    @test deₛ ≈ fd rtol=1e-6
    for Tc in 0.0:10.0:40.0
        T = Tc * u"°C"
        @test vapour_pressure(Bolton(), T) ≈ vapour_pressure(GoffGratch(), T) rtol=1e-2
    end
end
