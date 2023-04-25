using PyPlot
using Printf
using DelimitedFiles
fpath = pwd()
include(joinpath(fpath,"libs/BM_mod.jl"))


##
# params = Params(ϵ=0.005,ν=-1,φ=0.0*π/180,Da=0.0,w0=110.0,_hetero=true,dθ=1.38π/180)
params = Params(ϵ=0.001,ν=0.16,φ=50.0*π/180,Da=-0.0,w0=55.0,_hetero=true,dθ=0.7π/180)
initParamsWithStrain(params)
Latt = Lattice()
initLattice(Latt,params;lk=72)
blk = HBM()
initHBM(blk,Latt,params;lg=9)
writedlm(joinpath(fpath,"test/strain_bizhen.txt"),blk.Hk)
## 
# enegy jmap 
function plot_map(ϵ::Matrix{Float64},Latt::Lattice)
    fig,ax = subplots(2,figsize=(4,5))
    ϵ1 = reshape(ϵ[1,:],Latt.lk,Latt.lk)
    ϵ2 = reshape(ϵ[2,:],Latt.lk,Latt.lk)
    kvec = reshape(Latt.kvec,Latt.lk,Latt.lk) ./(params.kb)
    pl=ax[1].contourf(real(kvec),imag(kvec),ϵ1,20,cmap="coolwarm")
    colorbar(pl,ax=ax[1])
    ax[1].axis("equal")
    ax[1].plot([0;sqrt(3)/2],[1;1/2],"r+")

    pl=ax[2].contourf(real(kvec),imag(kvec),ϵ2,20,cmap="coolwarm_r")
    colorbar(pl,ax=ax[2])
    ax[2].axis("equal")
    ax[2].plot([0;sqrt(3)/2],[1;1/2],"r+")
    tight_layout()
    display(fig)
    # savefig(joinpath(fpath,"test/strained_BM_map.pdf"),transparent=true)
    close(fig)
end 

plot_map(readdlm(joinpath(fpath,"test/strain_bizhen.txt")),Latt)

## filling fraction map 
function plot_map_filling(ϵ::Matrix{Float64},Latt::Lattice)
    fig,ax = subplots(2,figsize=(4,5))
    ϵ1 = reshape(ϵ[1,:],Latt.lk,Latt.lk) 
    ϵ2 = reshape(ϵ[2,:],Latt.lk,Latt.lk)
    levels1 = collect(-0.9:0.1:-0.1) * maximum(abs.(ϵ))
    levels2 = collect(0.1:0.1:0.9) * maximum(abs.(ϵ))
    # levels2 = [0.247;0.328] * maximum(abs.(ϵ))
    ν1 = [8*sum( (sign.(levels1[i] .- ϵ[:] ) .+1 )./2 ) / length(ϵ[:]) - 4 for i in eachindex(levels1)] 
    ν2 = [8*sum( (sign.(levels2[i] .- ϵ[:] ) .+1)./2 ) / length(ϵ[:]) - 4 for i in eachindex(levels2)]
    
    # note that this definition of filling fraction is incorrect if maximum(ϵ1) > minimum(ϵ2)
    kvec = reshape(Latt.kvec,Latt.lk,Latt.lk) ./(params.kb)
    pl = ax[1].pcolor(real(kvec),imag(kvec),ϵ1,cmap="Blues_r")
    colorbar(pl,ax=ax[1])
    pl=ax[1].contour(real(kvec),imag(kvec),ϵ1,levels=levels1,cmap="tab10")
    ν1str = Dict(pl.levels[i]=> @sprintf("%1.1f",ν1[i]) for i in eachindex(pl.levels))
    ax[1].clabel(pl,pl.levels,fmt=ν1str,inline=true,fontsize=6)
    ax[1].axis("equal")
    ax[1].plot([0;sqrt(3)/2],[1;1/2],"r+")

    pl = ax[2].pcolor(real(kvec),imag(kvec),ϵ2,cmap="Greens")
    colorbar(pl,ax=ax[2])
    pl=ax[2].contour(real(kvec),imag(kvec),ϵ2,levels=levels2,cmap="tab10")
    ν2str = Dict(pl.levels[i]=> @sprintf("%1.1f",ν2[i]) for i in eachindex(pl.levels))
    
    ax[2].clabel(pl,pl.levels,fmt=ν2str,inline=true,fontsize=6)
    ax[2].axis("equal")
    ax[2].plot([0;sqrt(3)/2],[1;1/2],"r+")
    tight_layout()
    display(fig)
    close(fig)
end 

plot_map_filling(readdlm(joinpath(fpath,"test/strain_bizhen.txt")),Latt)


## filling fraction bounds van Hove 
function plot_map_filling_vanhove(ϵ::Matrix{Float64},Latt::Lattice)
    fig,ax = subplots(2,figsize=(4,5))
    ϵ1 = reshape(ϵ[1,:],Latt.lk,Latt.lk)
    ϵ2 = reshape(ϵ[2,:],Latt.lk,Latt.lk)
    # levels2 = [0.5435] * maximum(abs.(ϵ)) # 0.004, 40 degrees
    # levels2 = [0.4603,0.6219] * maximum(abs.(ϵ)) # 0.004, 30 degrees
    # levels2 = [0.0935,0.4475,0.555] * maximum(abs.(ϵ))  #0.003, 10 degrees
    levels1 = [-0.299,-0.281,-0.182] * maximum(abs.(ϵ)) # 0.000
    levels2 = [0.182,0.281,0.299] * maximum(abs.(ϵ)) # 0.000
    # levels2 = [0.3266,0.361,0.4286] * maximum(abs.(ϵ)) # 0.01, 0 degrees
    ν2 = [8*sum( (sign.(levels2[i] .- ϵ[:] ) .+1)./2 ) / length(ϵ[:]) - 4 for i in eachindex(levels2)]
    ν1 = [8*sum( (sign.(sort(levels1)[i] .- ϵ[:] ) .+1)./2 ) / length(ϵ[:]) - 4 for i in eachindex(levels2)]
    
    # note that this definition of filling fraction is incorrect if maximum(ϵ1) > minimum(ϵ2)
    kvec = reshape(Latt.kvec,Latt.lk,Latt.lk) ./(sqrt(3)*params.kb)

    pl = ax[1].contourf(real(kvec),imag(kvec),ϵ1,cmap="coolwarm",levels=20)
    colorbar(pl,ax=ax[1],ticks=-70:10:0)
    pl=ax[1].contour(real(kvec),imag(kvec),ϵ1,levels=sort(levels1),colors=["r","b","g"])
    ν2str = Dict(pl.levels[i]=> @sprintf("%1.2f",ν1[i]) for i in eachindex(pl.levels))
    ax[1].clabel(pl,pl.levels,fmt=ν2str,inline=true,fontsize=8)
    ax[1].plot(real([params.Kt+params.g1;params.Kb+params.g1+params.g2])./(sqrt(3)*params.kb),
            imag([params.Kt+params.g1;params.Kb+params.g1+params.g2])./(sqrt(3)*params.kb),"k*")
    ax[1].plot([0],[0],"kX")
    ax[1].axis("equal")

    pl = ax[2].contourf(real(kvec),imag(kvec),ϵ2,cmap="coolwarm_r",levels=20)
    colorbar(pl,ax=ax[2],ticks=0:10:70)
    pl=ax[2].contour(real(kvec),imag(kvec),ϵ2,levels=levels2,colors=["r","b","g"])
    ν2str = Dict(pl.levels[i]=> @sprintf("%1.2f",ν2[i]) for i in eachindex(pl.levels))
    ax[2].clabel(pl,pl.levels,fmt=ν2str,inline=true,fontsize=8)
    # ax[2].set_xticks(-0.5:0.5:1)
    ax[2].plot(real([params.Kt+params.g1;params.Kb+params.g1+params.g2])./(sqrt(3)*params.kb),
    imag([params.Kt+params.g1;params.Kb+params.g1+params.g2])./(sqrt(3)*params.kb),"k*")
    ax[2].plot([0],[0],"kX")
    ax[2].axis("equal")
    tight_layout()
    display(fig)
    close(fig)
end 

plot_map_filling_vanhove(readdlm(joinpath(fpath,"test/strain_bizhen.txt")),Latt)

##
# cut 
function plot_cuts(ϵ::Matrix{Float64},Latt::Lattice)
    lk = Latt.lk 
    ϵ = reshape(ϵ,:,lk,lk)
    # KΓ = [idK .- [2;1]*i for i in 0:(lk÷3)]
    # ΓM = [idΓ .+ [1;0]*i for i in 1:(lk÷2)]
    # MK = [idM .+ [1;2]*i for i in 1:(lk÷6)]
    idM = [lk÷2+1;1]
    idK = [2lk÷3+1;lk÷3+1]
    idK1 = [lk÷3+1;2lk÷3+1]
    idΓ = [1;1]

    fig = figure(figsize=(4,3))
    for iband in 1:size(ϵ,1)
        ϵK1K = [ϵ[iband,idK1[1]+i,idK1[2]-i] for i in 0:(lk÷3)]
        ϵKΓ = [ϵ[iband,idK[1]-2i,idK[2]-i] for i in 1:(lk÷3)]
        ϵΓM = [ϵ[iband,idΓ[1]+i,idΓ[2]] for i in 1:(lk÷2)]
        ϵcut = [ϵK1K;ϵKΓ;ϵΓM]
        plot(eachindex(ϵcut), ϵcut,"b-",ms=2)
    end
    xticks([1,(lk÷3+1),(lk÷3*2+1),(lk÷3*2+lk÷2+1)],
            [L"K_+",L"K_-",L"Γ",L"M"])
    # axvline(1)
    # axvline(length(ϵ1K1K))
    # axvline(length(ϵ1KΓ)+length(ϵ1K1K))
    # axvline(length(ϵ1KΓ)+length(ϵ1K1K)+length(ϵ1ΓM))
    title(L"w_0/w_1=1")
    # ylim([-6,6])
    tight_layout()
    display(fig)
    savefig(joinpath(fpath,"test/strain_ratio_10.pdf"),transparent=true)
    close(fig)
end
plot_cuts(blk.Hk,Latt)

