using PyPlot
using Printf
using DataFrames
using CSV
using DelimitedFiles
using JLD
# using Peaks
fpath = pwd()
include(joinpath(fpath,"libs/BM_mod.jl"))
include(joinpath(fpath,"libs/helpers_transportB0.jl"))

# basic settings
ϵ, ϕ, Da = 0.002, 0, -4100.0
str, str2 = Int(1000*ϵ), Int(Da)
lk = 256
params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da)
initParamsWithStrain(params)
Latt = Lattice()
initLattice(Latt,params;lk=lk)
blk = HBM()
# initHBM(blk,Latt,params;lg=9,_σrotation=false)
fname=joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).csv")
# df = DataFrame()
# df[!,"kx"] = [real(Latt.kvec[i]) for i in 1:length(Latt.kvec) for j in 1:2]
# df[!,"ky"] = [imag(Latt.kvec[i]) for i in 1:length(Latt.kvec) for j in 1:2]
# df[!,"Hk"] = blk.Hk[:]
# df[!,"Vx"] = blk.Vx[:]
# df[!,"Vy"] = blk.Vy[:]
# CSV.write(fname,df)
# df = DataFrame(CSV.File(fname))
blk.Hk = reshape(df[!,"Hk"],2,:)

# fname = joinpath(fpath,"strain/band_structure/eps_$(str)_Da_$(str2)_special_points_$(lk).csv")
# df = DataFrame()
# special_points, special_point_energies,special_point_fillings = 
#         get_special_points(ϵ=ϵ,ϕ=ϕ,Da=Da,_hetero=true,lk=lk,
#                         fname=joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).csv"))
# df[!,"$(ϕ)_k_R"] = real([special_points[1]; special_points[2]])
# df[!,"$(ϕ)_k_I"] = imag([special_points[1]; special_points[2]])
# df[!,"$(ϕ)_k_E"] = [special_point_energies[1]; special_point_energies[2]]
# df[!,"$(ϕ)_k_nu"] = [special_point_fillings[1]; special_point_fillings[2]]
# CSV.write(fname,df)

# compute transport coefficients
# transport_coefficients = mainTransportB0(blk,Latt,params,nμs=121);
# save(joinpath(fpath,"strain/transport/B0_eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).jld"),"transport",transport_coefficients)

transport_coefficients = load(joinpath(fpath,"strain/transport/B0_eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).jld"),"transport")


### Plots #### 
fname_special_points = joinpath(fpath,"strain/band_structure/eps_$(str)_Da_$(str2)_special_points_$(lk).csv")

# unsymmetrized transport
plot_ρs(transport_coefficients,fname_special_points)

# unsymmetrized derivative transport
plot_ρs_derivative(transport_coefficients,fname_special_points)

# symmetrization
νs,ρ,δθ = B0symmetrization(transport_coefficients) ;

# symmetrized transport
plot_symmetrized_transport(νs,ρ,fname_special_points)

# principal axis
plot_principle_strain_axis(νs,δθ,fname_special_points)
