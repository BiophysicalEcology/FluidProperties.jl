using FluidProperties
using Test
using Unitful

@testset "wet_bulb_properties" begin

    @testset "independent checks" begin
        P = 101325.0u"Pa"
        # Stull (2011) empirical formula, stated error within 1 K
        stull(T, RH) = let Tc = ustrip(u"°C", T), RH = 100RH
            Tc * atan(0.151977 * sqrt(RH + 8.313659)) + atan(Tc + RH) - atan(RH - 1.676331) +
                0.00391838 * RH^1.5 * atan(0.023101 * RH) - 4.686035
        end
        for Tc in 0.0:10.0:40.0, RH in 0.2:0.2:0.8
            T = Tc * u"°C"
            @test ustrip(u"°C", wet_bulb_temperature(T, RH, P)) ≈ stull(T, RH) atol=1.0
        end
        # The solution satisfies Davies-Jones eqn 2.3
        result = wet_bulb_properties(25.0u"°C", 0.5, P)
        x = (273.15u"K" / result.equivalent_temperature)^3.504
        @test FluidProperties._saturation_terms(result.wet_bulb_temperature, P).f ≈ x rtol=1e-5
        # FixedNewton agrees with RefinedNewton
        @test wet_bulb_temperature(30.0u"°C", 0.5, P; convergence=FixedNewton()) ≈
              wet_bulb_temperature(30.0u"°C", 0.5, P) atol=1e-3u"K"
        # Specific humidity equivalent to the relative humidity it implies
        rₛ = FluidProperties._saturation_terms(303.15u"K", P).rₛ
        @test wet_bulb_temperature(30.0u"°C", SpecificHumidity(0.5rₛ), P) ≈ wet_bulb_temperature(30.0u"°C", 0.5, P)
    end

    @testset "physical behaviour" begin
        P = 101325.0u"Pa"
        # Wet bulb is below air temperature and increases with humidity
        twb = [wet_bulb_temperature(25.0u"°C", rh, P) for rh in 0.1:0.1:0.9]
        @test issorted(twb)
        @test all(<(298.15u"K"), twb)
        # Kelvin and Celsius inputs agree
        @test wet_bulb_temperature(298.15u"K", 0.5, P) ≈ wet_bulb_temperature(25.0u"°C", 0.5, P)
        # Any compatible units are accepted
        @test wet_bulb_temperature(77.0u"°F", 0.5, 101.325u"kPa") ≈ wet_bulb_temperature(25.0u"°C", 0.5, P)
        @test wet_bulb_temperature(25.0u"°C", 0.5, 1013.25u"hPa") ≈ wet_bulb_temperature(25.0u"°C", 0.5, P)
        # Saturated air is at the air temperature
        @test wet_bulb_temperature(25.0u"°C", 1.0, P) ≈ 298.15u"K" atol=1e-3u"K"
    end

    @testset "units and types" begin
        result = wet_bulb_properties(25.0u"°C", 0.5, 101325.0u"Pa")
        @test result isa WetBulbProperties
        @test unit(result.wet_bulb_temperature) == u"K"
        @test unit(result.equivalent_temperature) == u"K"
        @test unit(result.equivalent_potential_temperature) == u"K"
        @test wet_bulb_temperature(25.0u"°C", 0.5, 101325.0u"Pa") == result.wet_bulb_temperature
        # Integer inputs
        @test wet_bulb_temperature(25u"°C", 1//2, 101325u"Pa") ≈ result.wet_bulb_temperature
    end

    @testset "invalid regime gives NaN" begin
        result = wet_bulb_properties(50.0u"°C", 0.9, 80000.0u"Pa")
        @test isnan(ustrip(result.wet_bulb_temperature))
        @test isnan(ustrip(result.equivalent_temperature))
        @test isfinite(ustrip(result.equivalent_potential_temperature))
    end

    @testset "humidity guards" begin
        P = 101325.0u"Pa"
        @test_throws DomainError wet_bulb_properties(25.0u"°C", 50.0, P)
        @test_throws DomainError wet_bulb_properties(25.0u"°C", -0.1, P)
        @test_throws DomainError wet_bulb_properties(25.0u"°C", SpecificHumidity(13.5), P)
        @test_throws DomainError wet_bulb_temperature(25.0u"°C", 50.0, P)
        @test isnan(ustrip(wet_bulb_temperature(25.0u"°C", NaN, P)))
    end

    @testset "missing input" begin
        @test ismissing(wet_bulb_properties(missing, 0.5, 101325.0u"Pa"))
        @test ismissing(wet_bulb_properties(25.0u"°C", missing, 101325.0u"Pa"))
        @test ismissing(wet_bulb_properties(25.0u"°C", 0.5, missing))
        @test ismissing(wet_bulb_temperature(missing, missing, missing))
    end

    @testset "type stable and non-allocating" begin
        T = 25.0u"°C"; P = 101325.0u"Pa"
        # Measured inside a function, so that the untyped loop variables below do not add allocations
        function allocations(T, humidity, P, convergence)
            wet_bulb_properties(T, humidity, P; convergence)
            return @allocated wet_bulb_properties(T, humidity, P; convergence)
        end
        for humidity in (0.5, SpecificHumidity(0.01)), convergence in (RefinedNewton(), FixedNewton())
            @inferred wet_bulb_properties(T, humidity, P; convergence)
            @inferred wet_bulb_temperature(T, humidity, P; convergence)
            @test allocations(T, humidity, P, convergence) == 0
        end
        @inferred wet_bulb_properties(50.0u"°C", 0.9, 80000.0u"Pa")
    end
end
