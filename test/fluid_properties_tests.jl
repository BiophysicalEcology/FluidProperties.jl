using FluidProperties
using Test
using Unitful

# Reference values are standard textbook properties (e.g. Incropera & DeWitt, Engineering
# Toolbox), with tolerances that allow for the regressions used by the package.

@testset "atmospheric_pressure" begin
    @test atmospheric_pressure(0.0u"m") ≈ 101325.0u"Pa"
    @test atmospheric_pressure(100.0u"m") == 100128.83927387102u"Pa"
    # International Standard Atmosphere at 5 km is 54 019 Pa
    @test atmospheric_pressure(5000.0u"m") ≈ 54019.0u"Pa" rtol=1e-3
    @test atmospheric_pressure(1000.0u"m"; reference_pressure=100000.0u"Pa") <
          atmospheric_pressure(1000.0u"m")
    @test ismissing(atmospheric_pressure(missing))
end

@testset "dry_air_properties" begin
    T = 20.0u"°C"
    P = 101325.0u"Pa"
    air = dry_air_properties(T, P)

    @test air.density ≈ 1.20u"kg/m^3" rtol=0.01
    @test air.dynamic_viscosity ≈ 1.81e-5u"kg/m/s" rtol=0.02
    @test air.kinematic_viscosity ≈ 1.51e-5u"m^2/s" rtol=0.02
    @test air.thermal_conductivity ≈ 0.0257u"W/m/K" rtol=0.02
    @test air.vapour_diffusivity ≈ 2.5e-5u"m^2/s" rtol=0.05
    @test uconvert(u"W/m^2", air.blackbody_emission) ≈ 418.7u"W/m^2" rtol=1e-3
    @test air.peak_wavelength ≈ 9.88e-6u"m" rtol=1e-3
    @test air.molar_mass ≈ 0.02885u"kg/mol" rtol=1e-3

    @testset "internal consistency" begin
        @test air.kinematic_viscosity ≈ air.dynamic_viscosity / air.density
        @test air.grashof_coefficient ≈ 9.80665u"m/s^2" / u"K"(T) / air.kinematic_viscosity^2 rtol=1e-6
        # Ideal gas law
        @test air.density ≈ air.molar_mass * P / (8.314462618u"J/mol/K" * u"K"(T)) rtol=1e-6
    end

    @testset "trends" begin
        hot = dry_air_properties(40.0u"°C", P)
        @test hot.density < air.density
        @test hot.dynamic_viscosity > air.dynamic_viscosity
        @test hot.thermal_conductivity > air.thermal_conductivity
        @test hot.peak_wavelength < air.peak_wavelength
        # Lower pressure lowers density and raises vapour diffusivity
        thin = dry_air_properties(T, 50000.0u"Pa")
        @test thin.density ≈ air.density * 50000 / 101325 rtol=1e-6
        @test thin.vapour_diffusivity > air.vapour_diffusivity
    end

    @testset "input variants" begin
        @test dry_air_properties(293.15u"K", P) == dry_air_properties(20.0u"°C", P)
        @test dry_air_properties(T; atmospheric_pressure=P) == air
        @test ismissing(dry_air_properties(missing, missing))
        @test ismissing(dry_air_properties(missing, P))
        @test ismissing(dry_air_properties(T, missing))
        @test ismissing(dry_air_properties(missing))
    end
end

@testset "wet_air_properties" begin
    T = 293.15u"K"
    P = 101325.0u"Pa"
    wet = wet_air_properties(T, 0.5, P)
    dry = wet_air_properties(T, 0.0, P)
    saturated = wet_air_properties(T, 1.0, P)
    P_sat = vapour_pressure(T)

    @test wet.vapour_pressure ≈ 0.5 * P_sat
    @test wet.relative_humidity == 0.5
    # Saturation vapour pressure at 20 °C is 2339 Pa
    @test saturated.vapour_pressure ≈ 2339.0u"Pa" rtol=1e-3
    # 50 % RH at 20 °C and sea level: mixing ratio ≈ 7.3 g/kg, vapour density ≈ 8.6 g/m³
    @test wet.mixing_ratio ≈ 0.0073 rtol=0.01
    @test wet.vapour_density ≈ 8.6e-3u"kg/m^3" rtol=0.01
    # Vapour density from the ideal gas law
    @test wet.vapour_density ≈ wet.vapour_pressure * 0.018015u"kg/mol" / (8.314462618u"J/mol/K" * T) rtol=0.01
    # Virtual temperature increment ≈ 0.61 r T
    @test wet.virtual_temp_increment ≈ 0.61 * wet.mixing_ratio * T rtol=0.05

    @testset "dry limit" begin
        @test dry.vapour_pressure == 0.0u"Pa"
        @test dry.mixing_ratio == 0.0
        @test dry.virtual_temp_increment == 0.0u"K"
        @test dry.specific_heat ≈ 1004.84u"J/kg/K"
        @test dry.water_potential == -999.0u"Pa"
        @test dry.density ≈ dry_air_properties(T, P).density rtol=2e-3
    end

    @testset "humidity trends" begin
        @test saturated.density < wet.density < dry.density
        @test dry.specific_heat < wet.specific_heat < saturated.specific_heat
        @test dry.mixing_ratio < wet.mixing_ratio < saturated.mixing_ratio
    end

    @testset "water potential" begin
        @test saturated.water_potential ≈ 0.0u"Pa" atol=1e-6u"Pa"
        # ψ = (R T / v_w) ln(RH) ≈ -93.8 MPa at 50 % RH and 20 °C
        @test wet.water_potential ≈ -93.8e6u"Pa" rtol=0.01
    end

    @testset "input variants" begin
        @test wet_air_properties(20.0u"°C", 0.5, P) == wet
        @test wet_air_properties(T, 0.5, P; vapour_pressure_equation=Huang()).vapour_pressure ≈
              0.5 * vapour_pressure(Huang(), T)
        @test ismissing(wet_air_properties(missing, missing, missing))
        @test ismissing(wet_air_properties(missing, 0.5, P))
        @test ismissing(wet_air_properties(T, missing, P))
        @test ismissing(wet_air_properties(T, 0.5, missing))
        @test ismissing(wet_air_properties(missing))
    end
end

@testset "water_properties" begin
    water = water_properties(20.0u"°C")
    # Handbook values for liquid water at 20 °C
    @test water.density ≈ 998.0u"kg/m^3" rtol=0.005
    @test water.specific_heat ≈ 4182.0u"J/kg/K" rtol=0.005
    @test water.thermal_conductivity ≈ 0.598u"W/m/K" rtol=0.01
    @test water.dynamic_viscosity ≈ 1.002e-3u"kg/m/s" rtol=0.05

    hot = water_properties(60.0u"°C")
    @test hot.density ≈ 983.2u"kg/m^3" rtol=0.005
    @test hot.thermal_conductivity ≈ 0.651u"W/m/K" rtol=0.01
    @test hot.dynamic_viscosity ≈ 4.66e-4u"kg/m/s" rtol=0.1
    @test hot.dynamic_viscosity < water.dynamic_viscosity
    @test hot.thermal_conductivity > water.thermal_conductivity

    # Kelvin input, and density is held at its 60 °C value above that
    @test water_properties(293.15u"K") == water
    @test water_properties(80.0u"°C").density == hot.density
    @test ismissing(water_properties(missing))
end

@testset "enthalpy_of_vaporisation" begin
    # Latent heat of vaporisation at 20 °C is 2454 kJ/kg
    @test enthalpy_of_vaporisation(20.0u"°C") ≈ 2454.0e3u"J/kg" rtol=1e-3
    @test enthalpy_of_vaporisation(40.0u"°C") ≈ 2406.0e3u"J/kg" rtol=2e-3
    @test enthalpy_of_vaporisation(293.15u"K") == enthalpy_of_vaporisation(20.0u"°C")
    # Sublimation below freezing is larger
    @test enthalpy_of_vaporisation(-10.0u"°C") > enthalpy_of_vaporisation(20.0u"°C")
    @test ismissing(enthalpy_of_vaporisation(missing))

    # Molar and mass-based values agree through the molar mass of water
    M_w = 0.018015u"kg/mol"
    @test molar_enthalpy_of_vaporisation(20.0u"°C") ≈ enthalpy_of_vaporisation(20.0u"°C") * M_w rtol=5e-3
    @test molar_enthalpy_of_vaporisation(293.15u"K") == molar_enthalpy_of_vaporisation(20.0u"°C")
    @test ismissing(molar_enthalpy_of_vaporisation(missing))
end

@testset "display of property structs" begin
    for properties in (
        dry_air_properties(20.0u"°C", 101325.0u"Pa"),
        wet_air_properties(20.0u"°C", 0.5, 101325.0u"Pa"),
        water_properties(20.0u"°C"),
        wet_bulb_properties(20.0u"°C", 0.5, 101325.0u"Pa"),
    )
        text = sprint(show, MIME("text/plain"), properties)
        @test startswith(text, string(nameof(typeof(properties))))
        @test !occursin("Quantity", text)
        @test count('
', text) == fieldcount(typeof(properties))
    end
    @test occursin("density", sprint(show, MIME("text/plain"), wet_air_properties(20.0u"°C", 0.5, 101325.0u"Pa")))
end
