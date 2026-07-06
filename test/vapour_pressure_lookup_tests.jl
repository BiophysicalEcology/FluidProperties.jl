using FluidProperties
import FluidProperties: VapourPressureLookup, GoffGratch, Teten, Huang
using Test
using Unitful

@testset "VapourPressureLookup" begin

    @testset "construction" begin
        # default construction succeeds
        vpl = VapourPressureLookup()
        @test vpl isa VapourPressureLookup

        # custom range and step
        vpl = VapourPressureLookup(; tmin=-80.0u"°C", tmax=100.0u"°C", step=0.05u"K")
        @test vpl isa VapourPressureLookup

        # alternative underlying formulations
        @test VapourPressureLookup(Teten()) isa VapourPressureLookup
        @test VapourPressureLookup(Huang()) isa VapourPressureLookup

        # step must be in Kelvin
        @test_throws ArgumentError VapourPressureLookup(step=0.1u"°C")
    end

    @testset "matches underlying formulation" begin
        # A fine step should make the lookup nearly identical to the
        # underlying equation at any in-range temperature.
        for formulation in (GoffGratch(), Huang(), Teten())
            vpl = VapourPressureLookup(formulation;
                tmin=-40.0u"°C", tmax=80.0u"°C", step=0.01u"K",
            )
            for Tc in -30.0:5.0:70.0
                T = (Tc + 273.15) * u"K"
                @test vapour_pressure(vpl, T) ≈ vapour_pressure(formulation, T) rtol=1e-4
            end
        end
    end

    @testset "returns Pa" begin
        vpl = VapourPressureLookup()
        @test unit(vapour_pressure(vpl, 293.15u"K")) == u"Pa"
    end

    @testset "linear interpolation at midpoint" begin
        # Halfway between two adjacent table nodes should equal the
        # average of the underlying equation evaluated at those nodes.
        formulation = GoffGratch()
        vpl = VapourPressureLookup(formulation;
            tmin=-10.0u"°C", tmax=40.0u"°C", step=1.0u"K",
        )
        # 290.15 K = 17°C and 291.15 K = 18°C land exactly on table nodes
        T1 = 290.15u"K"
        T2 = 291.15u"K"
        Tmid = 290.65u"K"
        es1 = vapour_pressure(formulation, T1)
        es2 = vapour_pressure(formulation, T2)
        @test vapour_pressure(vpl, Tmid) ≈ 0.5 * (es1 + es2) rtol=1e-10
    end

    @testset "exact at table endpoints" begin
        # At the tabulated nodes the lookup should return the underlying
        # equation exactly (no interpolation between adjacent cells).
        formulation = GoffGratch()
        vpl = VapourPressureLookup(formulation;
            tmin=0.0u"°C", tmax=30.0u"°C", step=1.0u"K",
        )
        for Tc in 1.0:5.0:25.0
            T = (Tc + 273.15) * u"K"
            @test vapour_pressure(vpl, T) ≈ vapour_pressure(formulation, T) rtol=1e-12
        end
    end

    @testset "missing input" begin
        vpl = VapourPressureLookup()
        @test ismissing(vapour_pressure(vpl, missing))
    end

end
