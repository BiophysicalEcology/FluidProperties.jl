using FluidProperties
using Test
using DataFrames
using CSV
using UnitfulMoles: @u_str

testdir = realpath(joinpath(dirname(pathof(FluidProperties)), "../test"))

# read in output from NicheMapR and input variables
pars_vec = DataFrame(CSV.File("$testdir/data/pars.csv"))[:, 2]
names = [:rh, :elevation]
# Zip into a NamedTuple
pars = (; zip(names, pars_vec)...)

T_airs = (DataFrame(CSV.File("$testdir/data/T_airs.csv"))[:, 2]) .* u"°C"
denair_dry = (DataFrame(CSV.File("$testdir/data/denair_dry.csv"))[:, 2]) .* u"kg/m^3"
visdyn = (DataFrame(CSV.File("$testdir/data/visdyn.csv"))[:, 2]) .* u"kg * m^-1 * s^-1"
viskin = (DataFrame(CSV.File("$testdir/data/viskin.csv"))[:, 2]) .* u"m^2 / s"
difvpr = (DataFrame(CSV.File("$testdir/data/difvpr.csv"))[:, 2]) .* u"m^2 / s"
thcond = (DataFrame(CSV.File("$testdir/data/thcond.csv"))[:, 2]) .* u"W * m^-1 * K^-1"
ggroup = (DataFrame(CSV.File("$testdir/data/ggroup.csv"))[:, 2]) .* u"m^-3 * K^-1"
bbemit = (DataFrame(CSV.File("$testdir/data/bbemit.csv"))[:, 2]) .* u"W/m^2"
emtmax = (DataFrame(CSV.File("$testdir/data/emtmax.csv"))[:, 2]) .* u"m"

esat = (DataFrame(CSV.File("$testdir/data/esat.csv"))[:, 2]) .* u"Pa"
e = (DataFrame(CSV.File("$testdir/data/e.csv"))[:, 2]) .* u"Pa"
vd = (DataFrame(CSV.File("$testdir/data/vd.csv"))[:, 2]) .* u"kg/m^3"
rw = (DataFrame(CSV.File("$testdir/data/rw.csv"))[:, 2])
tvinc = (DataFrame(CSV.File("$testdir/data/tvinc.csv"))[:, 2]) .* u"K"
denair_wet = (DataFrame(CSV.File("$testdir/data/denair_wet.csv"))[:, 2]) .* u"kg/m^3"
c_p_NMR = (DataFrame(CSV.File("$testdir/data/cp.csv"))[:, 2]) .* u"J/kg/K"
wtrpot = (DataFrame(CSV.File("$testdir/data/wtrpot.csv"))[:, 2]) .* u"Pa"

CP_water = (DataFrame(CSV.File("$testdir/data/CP_water.csv"))[:, 2]) .* u"J/kg/K"
DENSTY_water = (DataFrame(CSV.File("$testdir/data/DENSTY_water.csv"))[:, 2]) .* u"kg/m^3"
THCOND_water = (DataFrame(CSV.File("$testdir/data/THCOND_water.csv"))[:, 2]) .* u"W/m/K"
VISDYN_water = (DataFrame(CSV.File("$testdir/data/VISDYN_water.csv"))[:, 2]) .* u"kg * m^-1 * s^-1"

dry_air_out = dry_air_properties.(u"K".(T_airs), elevation=(pars.elevation)u"m")

# Extract each component into plain arrays
P_atmos           = getindex.(dry_air_out, 1)
ρ_air             = getindex.(dry_air_out, 2)
μ                 = getindex.(dry_air_out, 3)
ν                 = getindex.(dry_air_out, 4)
D_w               = getindex.(dry_air_out, 5)
k_air             = getindex.(dry_air_out, 6)
Grashof_group     = getindex.(dry_air_out, 7)
blackbody_emission = getindex.(dry_air_out, 8)
λ_max             = getindex.(dry_air_out, 9)

@testset "R DRYAIR comparisons" begin
    @test all(isapprox.(ρ_air, denair_dry; rtol=1e-9))
    @test all(isapprox.(μ, visdyn; rtol=1e-9))
    @test all(isapprox.(ν, viskin; rtol=1e-9))
    @test all(isapprox.(D_w, difvpr; rtol=1e-9))
    @test all(isapprox.(k_air, thcond; rtol=1e-9))
    @test all(isapprox.(Grashof_group, ggroup; rtol=1e-9))
    @test all(isapprox.(blackbody_emission, bbemit; rtol=1e-5))
    @test all(isapprox.(λ_max, emtmax; rtol=1e-9))
end 

# note FluidProperties.jl works with fractional relative humidity hence / 100.
wet_air_out = wet_air_properties.(u"K".(T_airs), P_atmos = P_atmos[1], rh = pars.rh / 100.)
# Extract each component into plain arrays
P_vap     = getindex.(wet_air_out, 1)
P_vap_sat = getindex.(wet_air_out, 2)
ρ_vap     = getindex.(wet_air_out, 3)
r_w       = getindex.(wet_air_out, 4)
T_vinc    = getindex.(wet_air_out, 5)
ρ_air_wet = getindex.(wet_air_out, 6)
c_p       = getindex.(wet_air_out, 7)
ψ         = getindex.(wet_air_out, 8)
rh        = getindex.(wet_air_out, 9)

@testset "R WETAIR comparisons" begin
    @test all(isapprox.(P_vap, e; rtol=1e-9))
    @test all(isapprox.(P_vap_sat, esat; rtol=1e-9))
    @test all(isapprox.(ρ_vap, vd; rtol=1e-9))
    @test all(isapprox.(r_w, rw; rtol=1e-9))
    @test all(isapprox.(T_vinc, tvinc; rtol=1e-9))
    @test all(isapprox.(ρ_air_wet, denair_wet; rtol=1e-9))
    @test all(isapprox.(c_p, c_p_NMR; rtol=1e-9))
    @test all(isapprox.(ψ, wtrpot; rtol=1e-9))
end 

water_properties_out = water_properties.(u"K".(T_airs))

# Extract each component into plain arrays
c_p_H2O           = getindex.(water_properties_out, 1)
ρ_H2O             = getindex.(water_properties_out, 2)
k_H2O             = getindex.(water_properties_out, 3)
μ_H2O             = getindex.(water_properties_out, 4)

@testset "R WATERPROP comparisons" begin
    @test all(isapprox.(c_p_H2O, CP_water; rtol=1e-9))
    @test all(isapprox.(ρ_H2O, DENSTY_water; rtol=1e-9))
    @test all(isapprox.(k_H2O, THCOND_water; rtol=1e-9))
    @test all(isapprox.(μ_H2O, VISDYN_water; rtol=1e-9))
end 