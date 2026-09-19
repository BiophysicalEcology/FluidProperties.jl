using FluidProperties
using Test
using Unitful

const P = 101325.0u"Pa"
const METHODS = (DaviesJones(), Stull(), Barenbrug(), Barenbrug(; vapour_pressure_equation=Huang()),
                 Smithsonian(), Smithsonian(; vapour_pressure_equation=Bolton()), EnergyBalance())
const ROOT_METHODS = (Barenbrug(), Smithsonian(), Smithsonian(; vapour_pressure_equation=Bolton()), EnergyBalance())

@testset "published values" begin
    # Stull (2011) example
    @test ustrip(u"°C", wet_bulb_temperature(Stull(), 20.0u"°C", 0.5, P)) ≈ 13.7 atol=0.05
    # List (1971) Table 98 example: 40.1 °C, 36.9 °C wet bulb and 923 hPa give a vapour pressure of 60.433 hPa
    RH_list = uconvert(NoUnits, 60.433u"hPa" / vapour_pressure(GoffGratch(), 40.1u"°C"))
    @test wet_bulb_temperature(Smithsonian(), 40.1u"°C", RH_list, 923.0u"hPa") ≈ u"K"(36.9u"°C") atol=0.02u"K"
    # Barenbrug (1948) example: 70 °F dry bulb, 60 °F wet bulb and 30 inHg give a vapour pressure of 0.4126 inHg
    inHg = 3386.389u"Pa"
    RH = 0.4126 / 0.7392  # saturation vapour pressure at 70 °F is 0.7392 inHg
    @test wet_bulb_temperature(Barenbrug(), 70.0u"°F", RH, 30 * inHg) ≈ u"K"(60.0u"°F") atol=0.05u"K"
end

@testset "agreement with Davies-Jones" begin
    for Tc in 0.0:10.0:50.0, RH in 0.2:0.2:0.8
        T = Tc * u"°C"
        reference = wet_bulb_temperature(DaviesJones(), T, RH, P)
        for method in ROOT_METHODS
            # Psychrometric equations differ from the pseudoadiabatic wet bulb by up to about 1 K at 50 °C
            @test wet_bulb_temperature(method, T, RH, P) ≈ reference atol=1.1u"K"
        end
        @test wet_bulb_temperature(Stull(), T, RH, P) ≈ reference atol=0.7u"K"
    end
end

@testset "psychrometric relations are satisfied" begin
    T = 30.0u"°C"
    RH = 0.4
    # Lemke and Kjellstrom (2012) eqn 4, with Bolton saturation vapour pressure
    T_w = wet_bulb_temperature(Smithsonian(; vapour_pressure_equation=Bolton()), T, RH, P)
    e = vapour_pressure(Bolton(), T_w) -
        0.00066 / u"K" * P * (u"K"(T) - T_w) * (1 + 0.00115 / u"K" * (T_w - 273.15u"K"))
    @test e ≈ RH * vapour_pressure(Bolton(), T) rtol=1e-6
    # Energy balance
    T_w = wet_bulb_temperature(EnergyBalance(), T, RH, P)
    @test T_w < u"K"(T)
    @test T_w > u"K"(T) - 30u"K"
end

@testset "limits" begin
    for method in ROOT_METHODS
        # Saturated air is at the air temperature
        @test wet_bulb_temperature(method, 25.0u"°C", 1.0, P) ≈ 298.15u"K" atol=1e-2u"K"
        # Wet bulb is below air temperature and rises with humidity
        twb = [wet_bulb_temperature(method, 25.0u"°C", RH, P) for RH in 0.0:0.25:0.75]
        @test issorted(twb)
        @test all(<(298.15u"K"), twb)
        # Lower pressure lowers the wet bulb
        @test wet_bulb_temperature(method, 25.0u"°C", 0.5, 50000.0u"Pa") <
              wet_bulb_temperature(method, 25.0u"°C", 0.5, P)
    end
    # Tighter tolerance and iteration limit
    @test wet_bulb_temperature(Barenbrug(; tolerance=1e-9u"K"), 25.0u"°C", 0.5, P) ≈
          wet_bulb_temperature(Barenbrug(), 25.0u"°C", 0.5, P) atol=1e-4u"K"
    @test isfinite(ustrip(wet_bulb_temperature(Barenbrug(; max_iterations=1), 25.0u"°C", 0.5, P)))
end

@testset "vapour pressure equation option" begin
    T = 25.0u"°C"
    @test wet_bulb_temperature(Smithsonian(; vapour_pressure_equation=Huang()), T, 0.5, P) ≈
          wet_bulb_temperature(Smithsonian(), T, 0.5, P) atol=0.05u"K"
    @test wet_bulb_temperature(EnergyBalance(; vapour_pressure_equation=Huang()), T, 0.5, P) ≈
          wet_bulb_temperature(EnergyBalance(), T, 0.5, P) atol=0.05u"K"
end

@testset "inputs" begin
    for method in METHODS
        twb = wet_bulb_temperature(method, 25.0u"°C", 0.5, P)
        @test unit(twb) == u"K"
        # Any compatible units
        @test wet_bulb_temperature(method, 77.0u"°F", 0.5, 101.325u"kPa") ≈ twb
        @test wet_bulb_temperature(method, 298.15u"K", 0.5, 1013.25u"hPa") ≈ twb
        # Humidity must be a fraction
        @test_throws DomainError wet_bulb_temperature(method, 25.0u"°C", 50.0, P)
        @test_throws DomainError wet_bulb_temperature(method, 25.0u"°C", -0.1, P)
        # Missing values
        @test ismissing(wet_bulb_temperature(method, missing, 0.5, P))
        @test ismissing(wet_bulb_temperature(method, 25.0u"°C", missing, P))
        @test ismissing(wet_bulb_temperature(method, 25.0u"°C", 0.5, missing))
    end
    # Specific humidity is only for Davies-Jones
    @test wet_bulb_temperature(DaviesJones(), 25.0u"°C", SpecificHumidity(0.01), P) isa Quantity
    @test_throws MethodError wet_bulb_temperature(Stull(), 25.0u"°C", SpecificHumidity(0.01), P)
    # Default method and keyword pass-through
    @test wet_bulb_temperature(25.0u"°C", 0.5, P) == wet_bulb_temperature(DaviesJones(), 25.0u"°C", 0.5, P)
    @test wet_bulb_temperature(25.0u"°C", 0.5, P; convergence=FixedNewton()) ==
          wet_bulb_temperature(DaviesJones(; convergence=FixedNewton()), 25.0u"°C", 0.5, P)
end

# Measured inside a function, so that the untyped loop variable does not add allocations
function allocations(method, T, RH, P)
    wet_bulb_temperature(method, T, RH, P)
    return @allocated wet_bulb_temperature(method, T, RH, P)
end

@testset "type stable and non-allocating" begin
    T = 25.0u"°C"
    for method in METHODS
        @inferred wet_bulb_temperature(method, T, 0.5, P)
        @test allocations(method, T, 0.5, P) == 0
    end
end

@testset "broadcasting" begin
    temperatures = [10.0, 20.0, 30.0] .* u"°C"
    for method in METHODS
        @test wet_bulb_temperature.(method, temperatures, 0.5, P) == [wet_bulb_temperature(method, T, 0.5, P) for T in temperatures]
    end
    @test wet_bulb_temperature.(DaviesJones(), temperatures, SpecificHumidity(0.01), P) isa Vector
end
