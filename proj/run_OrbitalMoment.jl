using PyPlot
using Printf
using DataFrames
using CSV
using DelimitedFiles
using JLD
fpath = pwd()
include(joinpath(fpath,"libs/BM_OrbitalMoment_mod.jl"))

params = Params(dθ=0.8π/180)
initParamsWithStrain(params)
Latt = Lattice()
initLattice(Latt,params;lk=64)
blk = HBM()
initHBM(blk,Latt,params;lg=9,_σrotation=false)
save("tmp_30_12.jld","blk",blk)
blk = load("tmp_30_08.jld","blk");
###
function plot_orbital_magnetic_moment(nband,blk,params)
    area_moire_k = imag(params.g1'*params.g2) 
    area_moire = imag(params.a1'*params.a2) 
    μs,νs,msr,mc,mtot = compute_orbital_magnetization(nband,blk,area_moire,area_moire_k)
    fig = figure(figsize=(4,3))
    plot(νs,-msr,"g--",label="SR")
    plot(νs,mc,"m--",label="Chern")
    plot(νs,mc-msr,"k-",label="Total")

    # data = readdlm("quantum_calculation.txt")
    # plot(data[:,1],data[:,2],"b-",label="Quantum")
    legend()
    axhline(0,ls=":",c="gray")
    xlabel("ν")
    ylabel(L"\rm M\ (μ_B\ per\ u.c.)")
    # ylim([-3,0.1])
    tight_layout()
    display(fig)
    # savefig("keep_30bands_64x64.pdf",transparent=true)
    close(fig)
    
    return nothing
end

plot_orbital_magnetic_moment(16,blk,params)

fig = figure()
ee = 1.602e-19 
hbar = 1.054571817e-34
aa = 2.46e-10
coeff = 1 
# pl = pcolor(reshape(real(Latt.kvec),Latt.lk,Latt.lk)/abs(params.g1),
#         reshape(imag(Latt.kvec),Latt.lk,Latt.lk)/abs(params.g1),
#        reshape(blk.Msrk[11,:],Latt.lk,Latt.lk)*ee*aa^2/hbar)
pl = pcolor(reshape(real(Latt.kvec),Latt.lk,Latt.lk)/abs(params.g1),
        reshape(imag(Latt.kvec),Latt.lk,Latt.lk)/abs(params.g1),
        log.(abs.(reshape(blk.Ωnk[16,:]*coeff,Latt.lk,Latt.lk))))
colorbar(pl)
axis("equal")
tight_layout()
display(fig)
close(fig)

sum(blk.Ωnk[16,:]) * imag(params.g1'*params.g2)  / length(blk.Ωnk[16,:]) /(2π) 