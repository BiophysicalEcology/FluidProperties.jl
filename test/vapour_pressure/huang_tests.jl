using FluidProperties
import FluidProperties: Huang
using Test
using Unitful

# Reference values from Tables 1 and 2 in Huang (2018),
# "A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice",
# J. Appl. Meteor. Climatol. 57, 1265–1272.
# These are IAPWS reference values; the Huang formula has MAXRE 0.0057% (water) and 0.023% (ice).

@testset "Huang (2018) against IAPWS reference values" begin

    # Table 1: saturation vapour pressure over liquid water (t > 0°C)
    @testset "liquid water" begin
        water_cases = [
            (0.01,  611.655),
            (20.0,  2339.32),
            (40.0,  7384.94),
            (60.0,  19946.4),
            (80.0,  47414.5),
            (100.0, 101418.0),
        ]
        for (t_C, es_ref) in water_cases
            T = (t_C + 273.15) * u"K"
            @test vapour_pressure(Huang(), T) ≈ es_ref * u"Pa" rtol=0.0001
        end
    end

    # Table 2: saturation vapour pressure over ice (t ≤ 0°C)
    @testset "ice" begin
        ice_cases = [
            (-100.0, 0.0014049),
            (-80.0,  0.054773),
            (-60.0,  1.0813),
            (-40.0,  12.8412),
            (-20.0,  103.239),
            (0.0,    611.153),
        ]
        for (t_C, es_ref) in ice_cases
            T = (t_C + 273.15) * u"K"
            @test vapour_pressure(Huang(), T) ≈ es_ref * u"Pa" rtol=0.0003
        end
    end

end
