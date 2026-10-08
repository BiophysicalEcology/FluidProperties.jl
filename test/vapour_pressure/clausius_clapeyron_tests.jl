using FluidProperties
using Test
using Unitful

# Reference values from Thermodynamics.jl v1.3.0 with the default ClimaParams.jl parameters:
# (T in °C, saturation_vapor_pressure(param_set, T, Liquid()), saturation_vapor_pressure(param_set, T, Ice()))
const THERMODYNAMICS_REFERENCE = [
    (-40.0, 19.0171017998441, 12.829387821882092),
    (-20.0, 125.61075295754343, 103.29543988822948),
    (-10.0, 286.57899125208445, 259.9734334778779),
    (0.0, 611.2129404468477, 611.1537306200953),
    (0.01, 611.6569999999989, 611.6569999999984),
    (10.0, 1227.6651908618685, 1351.7436550814323),
    (20.0, 2337.0605652900877, 2830.6031357330617),
    (30.0, 4239.89436458799, 5642.395558133531),
    (40.0, 7365.740643349435, 10757.349043575327),
    (50.0, 12304.753320247039, 19697.284880513573),
]

@testset "Clausius-Clapeyron vapour pressure" begin
    liquid = ClausiusClapeyron(; ice=false)
    @testset "matches Thermodynamics.jl" begin
        for (Tc, e_liquid, e_ice) in THERMODYNAMICS_REFERENCE
            T = Tc * u"°C"
            @test vapour_pressure(liquid, T) ≈ e_liquid * u"Pa" rtol=1e-12
            # Ice below the triple point, liquid above
            expected = u"K"(T) < 273.16u"K" ? e_ice : e_liquid
            @test vapour_pressure(ClausiusClapeyron(), T) ≈ expected * u"Pa" rtol=1e-12
        end
    end
    @test vapour_pressure(ClausiusClapeyron(), 273.16u"K") ≈ 611.657u"Pa"
    @test unit(vapour_pressure(ClausiusClapeyron(), 20.0u"°C")) == u"Pa"
    @test ismissing(vapour_pressure(ClausiusClapeyron(), missing))
    for Tc in -40.0:10.0:50.0
        T = Tc * u"°C"
        @test vapour_pressure(ClausiusClapeyron(), T) ≈ vapour_pressure(GoffGratch(), T) rtol=1e-2
    end
    @test (@allocated vapour_pressure(ClausiusClapeyron(), 293.15u"K")) == 0
end
