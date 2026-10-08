using FluidProperties
using Test
using Unitful

@testset "natural_wet_bulb_temperature" begin
    T = 30.0u"°C"
    T_w = 22.0u"°C"
    ΔT = 8u"K"
    T_K = u"K"(T)

    @testset "no radiant heat" begin
        # C = 0.85 when still, 1 when ventilated, logarithmic between
        @test natural_wet_bulb_temperature(T, T_w, 0.0u"m/s") ≈ T_K - 0.85 * ΔT
        @test natural_wet_bulb_temperature(T, T_w, 0.02u"m/s") ≈ T_K - 0.85 * ΔT
        @test natural_wet_bulb_temperature(T, T_w, 1.0u"m/s") ≈ T_K - 0.96 * ΔT
        @test natural_wet_bulb_temperature(T, T_w, 2.0u"m/s") ≈ T_K - (0.96 + 0.069 * log10(2)) * ΔT
        @test natural_wet_bulb_temperature(T, T_w, 4.0u"m/s") == u"K"(T_w)
        # Continuous at the limits, within the rounding of the fit
        @test natural_wet_bulb_temperature(T, T_w, 0.031u"m/s") ≈ natural_wet_bulb_temperature(T, T_w, 0.029u"m/s") atol=0.1u"K"
        @test natural_wet_bulb_temperature(T, T_w, 3.01u"m/s") ≈ natural_wet_bulb_temperature(T, T_w, 2.99u"m/s") atol=0.1u"K"
        # Natural wet bulb is above the psychrometric value
        @test natural_wet_bulb_temperature(T, T_w, 0.5u"m/s") > u"K"(T_w)
        # Globe less than 4 K above air is ignored
        @test natural_wet_bulb_temperature(T, T_w, 4.0u"m/s"; globe_temperature=32.0u"°C") == u"K"(T_w)
    end

    @testset "radiant heat" begin
        Tg = 40.0u"°C"
        @test natural_wet_bulb_temperature(T, T_w, 0.05u"m/s"; globe_temperature=Tg) ≈ u"K"(T_w) + 2.5u"K" + 1.1u"K"
        @test natural_wet_bulb_temperature(T, T_w, 2.0u"m/s"; globe_temperature=Tg) ≈ u"K"(T_w) + 2.5u"K" - 0.1u"K"
        @test natural_wet_bulb_temperature(T, T_w, 0.5u"m/s"; globe_temperature=Tg) ≈
              u"K"(T_w) + 2.5u"K" + 0.1u"K" / 0.5^1.1 - 0.2u"K"
        @test natural_wet_bulb_temperature(T, T_w, 1.0u"m/s"; globe_temperature=Tg) ≈ u"K"(T_w) + 2.5u"K" - 0.1u"K"
        # Threshold is 4 K
        @test natural_wet_bulb_temperature(T, T_w, 2.0u"m/s"; globe_temperature=34.0u"°C") ≈
              u"K"(T_w) + 1.0u"K" - 0.1u"K"
    end

    @test unit(natural_wet_bulb_temperature(T, T_w, 1.0u"m/s")) == u"K"
    @test natural_wet_bulb_temperature(T, T_w, 3.6u"km/hr") ≈ natural_wet_bulb_temperature(T, T_w, 1.0u"m/s")
    @test_throws DomainError natural_wet_bulb_temperature(T, T_w, -1.0u"m/s")
    @test ismissing(natural_wet_bulb_temperature(missing, T_w, 1.0u"m/s"))
    @test ismissing(natural_wet_bulb_temperature(T, T_w, missing; globe_temperature=T))

    @inferred natural_wet_bulb_temperature(T, T_w, 1.0u"m/s")
    @inferred natural_wet_bulb_temperature(T, T_w, 1.0u"m/s"; globe_temperature=40.0u"°C")
    function allocations(T, T_w, v; kw...)
        natural_wet_bulb_temperature(T, T_w, v; kw...)
        return @allocated natural_wet_bulb_temperature(T, T_w, v; kw...)
    end
    @test allocations(T, T_w, 1.0u"m/s") == 0
    @test allocations(T, T_w, 1.0u"m/s"; globe_temperature=T + 10u"K") == 0
end
