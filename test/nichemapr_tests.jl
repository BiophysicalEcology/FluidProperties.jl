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

P_atmos = atmospheric_pressure((pars.elevation)u"m")
dry_air_out = dry_air_properties.(u"K".(T_airs), P_atmos)

# Extract each component into plain arrays
ρ_air             = getfield.(dry_air_out, :ρ_air)
μ                 = getfield.(dry_air_out, :μ)
ν                 = getfield.(dry_air_out, :ν)
D_w               = getfield.(dry_air_out, :D_w)
k_air             = getfield.(dry_air_out, :k_air)
Grashof_group     = getfield.(dry_air_out, :Grashof_group)
blackbody_emission = getfield.(dry_air_out, :blackbody_emission)
λ_max             = getfield.(dry_air_out, :λ_max)

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
wet_air_out = wet_air_properties.(u"K".(T_airs), pars.rh / 100., P_atmos)

# Extract each component into plain arrays
P_vap     = getfield.(wet_air_out, :P_vap)
ρ_vap     = getfield.(wet_air_out, :ρ_vap)
r_w       = getfield.(wet_air_out, :r_w)
T_vinc    = getfield.(wet_air_out, :T_vinc)
c_p       = getfield.(wet_air_out, :c_p)
ψ         = getfield.(wet_air_out, :ψ)
rh        = getfield.(wet_air_out, :rh)

@testset "R WETAIR comparisons" begin
    @test all(isapprox.(P_vap, e; rtol=1e-9))
    @test all(isapprox.(ρ_vap, vd; rtol=1e-9))
    @test all(isapprox.(r_w, rw; rtol=1e-9))
    @test all(isapprox.(T_vinc, tvinc; rtol=1e-9))
    @test all(isapprox.(c_p, c_p_NMR; rtol=1e-9))
    @test all(isapprox.(ψ, wtrpot; rtol=1e-9))
end 

water_properties_out = water_properties.(u"K".(T_airs))

# Extract each component into plain arrays
c_p_H2O           = getindex.(water_properties_out, :c_p_H2O)
ρ_H2O             = getindex.(water_properties_out, :ρ_H2O)
k_H2O             = getindex.(water_properties_out, :k_H2O)
μ_H2O             = getindex.(water_properties_out, :μ_H2O)

@testset "R WATERPROP comparisons" begin
    @test all(isapprox.(c_p_H2O, CP_water; rtol=1e-9))
    @test all(isapprox.(ρ_H2O, DENSTY_water; rtol=1e-9))
    @test all(isapprox.(k_H2O, THCOND_water; rtol=1e-9))
    @test all(isapprox.(μ_H2O, VISDYN_water; rtol=1e-9))
end 