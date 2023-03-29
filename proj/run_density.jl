using PyPlot
using Printf
using DataFrames
using CSV
using DelimitedFiles
fpath = pwd()
include(joinpath(fpath,"libs/BM_mod.jl"))

##
lk = 64
ϵs = collect(0.000:0.001:0.007)
ϕs = collect(0:2:60)
Da= -4100.0
# Da = 0.0
iϵ =  1 #parse(Int,ARGS[1])
iϕ =  1 #parse(Int,ARGS[2])
ϵ = 0.01 #ϵs[iϵ]
ϕ = ϕs[iϕ]
params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da,_hetero=true,dθ=1.38π/180)
initParamsWithStrain(params)
Latt = Lattice()
initLattice(Latt,params;lk=lk)
blk = HBM()
initHBM(blk,Latt,params;lg=9,_σrotation=false)


Uk = reshape(blk.Uk,4,blk.lg^2,2,length(Latt.kvec))

rvec = reshape(range(-2,2,61),:,1)*params.a1 .+ reshape(range(-2,2,61),1,:) * params.a2
ρvec = zeros(4,size(rvec,1),size(rvec,2))

for y in 1:size(rvec,2), x in 1:size(rvec,1)
    r = rvec[x,y]
    expigr = reshape(exp.(1im*real(r'.*blk.gvec)),1,:,1,1)
    expikr = reshape(exp.(1im*real(r'.*Latt.kvec)),1,1,1,:)
    ψnkr =  reshape(sum( expigr .* expikr .* Uk, dims=(2)),(4,2,:)) / sqrt(length(Latt.kvec)*length(blk.gvec))
    ρvec[:,x,y] = reshape(sum(abs2.(ψnkr),dims=(2,3)),:)
end

fig = figure()
pcolormesh(real(rvec),imag(rvec),(reshape(sum(ρvec,dims=1),size(rvec,1),:)))
# for r in -2:2
#     axhline(r)
# end
# for c in -2:2 
#     axvline(c)
# end
arrow(0,0,real(params.a1),imag(params.a1))
arrow(0,0,real(params.a2),imag(params.a2))
arrow(0,0,real(params.a2-params.a1),imag(params.a2-params.a1))
xlabel(L"x/|L_1|")
ylabel(L"y/|L_1|")
# title("ϵ=1%,φ=0")
title("unstrained")
axis("equal")
tight_layout()
display(fig)
savefig("dos_unstrained.pdf",transparent=true)
close(fig)



kvec = reshape(Latt.kvec,64,64)
fig = figure()
contourf(real(kvec),imag(kvec),reshape(blk.Hk[2,:],64,64))
# for r in -2:2
#     axhline(r)
# end
# for c in -2:2 
#     axvline(c)
# end
axis("equal")
tight_layout()
display(fig)
savefig("energy_strained.pdf",transparent=true)
close(fig)