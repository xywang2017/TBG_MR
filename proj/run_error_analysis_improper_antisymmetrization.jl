using PyPlot
using Printf
using DataFrames
using CSV
using DelimitedFiles
using JLD
# using Peaks
fpath = pwd()
include(joinpath(fpath,"libs/BM_mod.jl"))
include(joinpath(fpath,"libs/helpers_ASM.jl"))
include(joinpath(fpath,"libs/helpers.jl"))
# basic settings
ϵ, ϕ, Da = 0.003, 0, -4100.0
str, str2 = Int(1000*ϵ), Int(Da)
lk = 384
params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da)
initParamsWithStrain(params)
Latt = Lattice()
initLattice(Latt,params;lk=lk)
blk = HBM()
fname=joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).csv")
df = DataFrame(CSV.File(fname))
blk.Hk = reshape(df[!,"Hk"],2,:)
B0 = 4
eBτs = B0*[1;-1;-1+0.025;-1-0.025;-1+0.05;-1-0.05;-1+0.1;-1-0.1;-1+0.2;-1-0.2;-1+0.4;-1-0.4].*1e-4
# transport_coefficients = mainTransport(blk,Latt,params;eBτs = eBτs,nμs=121);
# save(joinpath(fpath,"strain/transport/ASM_B0$(B0)_eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).jld"),"transport",transport_coefficients)

ϵa = 1260.88
ω0τs = eBτs * ϵa
transport_coefficients = load(joinpath(fpath,"strain/transport/ASM_B0$(B0)_eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).jld"),"transport")

### Plots #### 
fname_special_points = joinpath(fpath,"strain/band_structure/eps_$(str)_Da_$(str2)_special_points.csv")

# symmetrization
νs,ρ,δρ0 = Bsymmetrization(transport_coefficients,eBτs,0,[1,2]) ;
νs,ρ,δρ1 = Bsymmetrization(transport_coefficients,eBτs,0*π/180,[1,4]) ;
νs,ρ,δρ2 = Bsymmetrization(transport_coefficients,eBτs,0*π/180,[1,7]) ;

plot_hall(νs,[δρ0,δρ2],fname_special_points)