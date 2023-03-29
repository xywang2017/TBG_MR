using PyPlot
using Printf
using DataFrames
using CSV
using DelimitedFiles
fpath = pwd()
include(joinpath(fpath,"libs/BM_mod.jl"))
include(joinpath(fpath,"libs/helpers.jl"))

##
lk = 128
ϵs = collect(0.001:0.001:0.007)
ϕs = collect(0:2:60)
Da= -4100.0
iϵ =  2 #parse(Int,ARGS[1])
iϕ =  1 #parse(Int,ARGS[2])
ϵ = ϵs[iϵ]
ϕ = ϕs[iϕ]
# println(ϵ," ",ϕ," start")
str = Int(1000 * ϵ)
str2 = Int(Da)
# params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da,_hetero=true,dθ=1.38π/180)
params = Params(ϵ=0.0,φ=ϕ*π/180,Da=Da,_hetero=true,dθ=1.39π/180)
initParamsWithStrain(params)
# get_bandstructure_info(params;ϵ=ϵ,ϕ=ϕ,Da=Da,_hetero=true,lk=lk,
#         fname=joinpath(fpath,"strain/band_structure/w0_0_eps_$(str)_phi_$(ϕ)_Da_$(str2).csv"))
# println(ϵ," ",ϕ," end")
# ----

fname = joinpath(fpath,"strain/band_structure/eps_$(str)_Da_$(str2)_special_points.csv")
df = DataFrame()
for ϕ in ϕs
    special_points, special_point_energies,special_point_fillings = 
        get_special_points(params;ϵ=ϵ,ϕ=ϕ,Da=Da,_hetero=true,lk=lk,
                        fname=joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2).csv"))
    if length(special_points[1]) + length(special_points[2]) != 12
        println("number of special point are: ",length(special_points[1])," ",length(special_points[2]))
        ntotal = length(special_points[1])+length(special_points[2])
        df[!,"$(ϕ)_k_R"] = real([special_points[1]; special_points[2]])
        df[!,"$(ϕ)_k_I"] = imag([special_points[1]; special_points[2]])
        df[!,"$(ϕ)_k_E"] = [special_point_energies[1]; special_point_energies[2]]
        df[!,"$(ϕ)_k_nu"] = [special_point_fillings[1]; special_point_fillings[2]]
    else
        df[!,"$(ϕ)_k_R"] = real([special_points[1]; special_points[2]])
        df[!,"$(ϕ)_k_I"] = imag([special_points[1]; special_points[2]])
        df[!,"$(ϕ)_k_E"] = [special_point_energies[1]; special_point_energies[2]]
        df[!,"$(ϕ)_k_nu"] = [special_point_fillings[1]; special_point_fillings[2]]
    end
end
CSV.write(fname,df)

# --------------------- example of tilted Dirac cone --------------- # 
# function tilted_dirac_cone_ex()
#     Da= -4100.0
#     ϵ = 0.005
#     ϕ = 40
#     str = Int(1000 * ϵ)
#     str2 = Int(Da)
#     params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da,_hetero=true)
#     initParamsWithStrain(params)
#     fname = joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2).csv")
#     Latt = Lattice()
#     initLattice(Latt,params;lk=128)
#     df = DataFrame(CSV.File(fname))
#     energies = reshape(df[!,"Hk"],2,128,128)
#     kvec = reshape(df[!,"kx"]+1im*df[!,"ky"],2,128,128)[1,:,:]

#     special_points, special_point_energies,special_point_fillings = 
#         get_special_points(ϵ=ϵ,ϕ=ϕ,Da=Da,_hetero=true,lk=128,
#                         fname=joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2).csv"))
#     contours = [31.528508442825924,31.941471278432555]
#     contour_colors=["b","r"]
#     special_points = special_points[2][sortperm(special_point_energies[2])]

#     fig = figure(figsize=(3,2))
#     νs = [2.29638671875,2.3798828125]
#     pl = contourf(real(kvec),imag(kvec),energies[2,:,:],levels=20,cmap="Spectral_r",vmin=-10,vmax=70)
#     # colorbar(pl,ticks=collect(0:10:60))

#     if !isempty(contours)
#         pl = contour(real(kvec),imag(kvec),energies[2,:,:],levels=contours,colors=contour_colors,linestyles="solid")
#         contour_labels = Dict(pl.levels[i]=> @sprintf("%1.1f",νs[i]) for i in eachindex(pl.levels))
#         clabel(pl,pl.levels,fmt=contour_labels,inline=true,fontsize=6)
#     end
#     DP1 = 0.11159472332552663 + 0.012430432107548768im
#     DP2 = 0.022756270233947103 + 0.10338800844860978im
#     plot(real(special_points[1]),imag(special_points[1]),"kX")
#     plot(real(special_points[2]),imag(special_points[2]),"o",c="b")
#     plot(real(special_points[3]),imag(special_points[3]),"o",c="r")
#     plot(real(special_points[4]),imag(special_points[4]),"kX")
#     plot(real(DP1),imag(DP1),"k*")
#     plot(real(DP2),imag(DP2),"k*")
#     axis("off")
#     axis("equal")
#     tight_layout()
#     savefig("contours.pdf",transparent="true")
#     display(fig)
#     close(fig)


    
#     ## plot map and cuts over the Dirac point for a given strain ### 
#     dp = 0.022756270233947103 + 0.10338800844860978im
#     dp1, dp2 = real(dp'*params.a1)/(2π) , real(dp'*params.a2)/(2π)
#     kvec = collect(range(0,1,length=60))*params.g1 .+ dp2 *params.g2
#     latt = initLatticeWithKvec(kvec)
#     blk = HBM()
#     initHBM(blk,latt,params;lg=9)
#     fig = figure(figsize=(3,2.5))
#     dk = Float64[0]
#     for i in 2:length(kvec)
#         push!(dk,dk[i-1]+abs(kvec[i]-kvec[i-1]))
#     end
#     plot(dk/abs(params.g1),blk.Hk[1,:],"r-",label="lower band",lw=2)
#     plot(dk/abs(params.g1),blk.Hk[2,:],"b-",label="upper band",lw=2)
#     legend()
#     # plot(eachindex(data),data[sortperm(data)],"bo")
#     ylabel("E (meV)")
#     xlabel(L"k_1")
#     tight_layout()
#     savefig("tilted_DP.pdf",transparent=true)
#     display(fig)
#     close(fig)
# end

# --------------------- plot special points --------------- # 
# fname = joinpath(fpath,"strain/band_structure/eps_3_Da_-4100_special_points.csv")
# savename=joinpath(fpath,"strain/band_structure/eps_3_Da_-4100_special_points_filling.pdf")
# plot_special_points(fname=fname,flag="filling",savename=savename)
# plot_special_points_combined(fname=fname,flag="energy",savename=savename)

## plot bounds 
# function plot_bounds()
#     df = DataFrame(CSV.File(joinpath(fpath,"alpha0.8/bounds.csv")))
#     fig, ax = subplots(1,2,figsize=(5,2.8))
#     ax[1].plot(df[!,"epsilon"],df[!,"delta_max"],"b<-",label=L"Δ_{max}")
#     ax[1].plot(df[!,"epsilon"],df[!,"delta_min"],"r^-",label=L"Δ_{min}")
#     ax[1].set_xlabel("ϵ (%)")
#     ax[1].set_ylabel("Δ (meV)")
#     ax[1].legend()

#     ax[2].plot(df[!,"epsilon"],df[!,"dnu_upper"],"b<-",label="upper")
#     ax[2].plot(df[!,"epsilon"],df[!,"dnu_lower"],"r^-",label="lower")
#     ax[2].set_xlabel("ϵ (%)")
#     ax[2].set_ylabel("δν (open FS)")
#     ax[2].legend()
    
#     tight_layout()
#     savefig("bounds.pdf",transparent=true)
#     display(fig)
#     close(fig)
# end
# plot_bounds()

# --------------------- plot energy maps --------------- # 
lk = 128
ϵ = 0.002
str = Int(1000 * ϵ)
for ϕ in 0:0
    params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da,_hetero=true)
    initParamsWithStrain(params)
    Latt = Lattice()
    initLattice(Latt,params;lk=lk)
    blk = HBM()
    df = DataFrame(CSV.File(joinpath(fpath,"strain/band_structure/w0_110_eps_$(str)_phi_$(ϕ)_Da_$(str2).csv")))
    blk.Hk = reshape(df[!,"Hk"],2,:)
    kvec = reshape(Latt.kvec,Latt.lk,Latt.lk) ./ (sqrt(3)*params.kb)
    df = DataFrame( CSV.File(joinpath(fpath,"strain/band_structure/eps_$(str)_Da_$(str2)_special_points.csv")) ) 
    data = df[!,"$(ϕ)_k_E"]
    special_points = df[!,"$(ϕ)_k_R"] + 1im * df[!,"$(ϕ)_k_I"]
    special_points = special_points[sortperm(data)] ./(sqrt(3)*params.kb)
    data = data[sortperm(data)]
    νs = [4*sum((sign.(data[i].-blk.Hk).+1)./2)/Latt.lk^2 - 4 for i in 9:11]
    savename = joinpath(fpath,"strain/band_structure/map_eps_$(str)_phi_$(ϕ)_Da_$(str2).pdf")
    plot_contour_maps(real(kvec),imag(kvec),reshape(blk.Hk,2,Latt.lk,Latt.lk)[2,:,:];
                contours=data[9:11],contour_colors=["b","r","g"],νs=νs,
                special_points=special_points[[6;8;9:12]], figsize=(3,2.), _legend=false, cmap="Spectral_r",
                savename = savename)
    # νs = [4*sum((sign.(data[i].-blk.Hk).+1)./2)/Latt.lk^2 - 4 for i in 9:9]
    # plot_contour_maps(real(kvec),imag(kvec),reshape(blk.Hk,2,Latt.lk,Latt.lk)[2,:,:];
    #             contours=data[9:9],contour_colors=["gray"],νs=νs,
    #             special_points=special_points[[6;8;9:11;14]], figsize=(3,2.), _legend=false, cmap="Spectral_r",
    #             savename = savename)
end

# --------------------- plot cut --------------- # 
lk = 128
for ϕ in 0:0
    params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da,_hetero=true)
    initParamsWithStrain(params)
    Latt = Lattice()
    initLattice(Latt,params;lk=lk)
    blk = HBM()
    # df = DataFrame(CSV.File(joinpath(fpath,"strain/band_structure/w0_0_eps_$(str)_phi_$(ϕ)_Da_$(str2).csv")))
    df = DataFrame(CSV.File(joinpath(fpath,"strain/band_structure/w0_0_eps_2_phi_$(ϕ)_Da_$(str2).csv")))
    blk.Hk = reshape(df[!,"Hk"],2,:)
    energies = reshape(blk.Hk,2,Latt.lk,Latt.lk)
    kvec = reshape(Latt.kvec,Latt.lk,Latt.lk) ./ (sqrt(3)*params.kb)

    fig = figure(figsize=(4,3))
    plot(real(kvec[:,1]),energies[2,:,1],"m-",label="cut g1")
    plot(real(kvec[:,1]),energies[2,1,:],"g-",label="cut g2")
    plot(real(kvec[:,1]),diag(energies[2,:,:]),"b-",label="cut g2+g2")
    xlabel("k")
    ylabel("E")
    legend()
    tight_layout()
    display(fig)
    close(fig)
end

# --------------------- plot velocities --------------- # 
lk = 128
for ϕ in 40:40
    params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da,_hetero=true)
    initParamsWithStrain(params)
    Latt = Lattice()
    initLattice(Latt,params;lk=lk)
    blk = HBM()
    df = DataFrame(CSV.File(joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2).csv")))
    blk.Hk = reshape(df[!,"Hk"],2,:)
    blk.Vx = reshape(df[!,"Vx"],2,:)
    blk.Vy = reshape(df[!,"Vy"],2,:)
    velocity = reshape(sqrt.(blk.Vx.^2+blk.Vy.^2), 2,Latt.lk,Latt.lk)
    kvec = reshape(Latt.kvec,Latt.lk,Latt.lk) ./ (sqrt(3)*params.kb)

    fig = figure(figsize=(4,3))
    kvals = (real(kvec[2:end,1]) + real(kvec[1:(end-1),1]) ) /2
    dk = diff(real(kvec[:,1]))
    plot(kvals,diff(velocity[2,:,1])./dk,"m-",label="cut g1")
    plot(kvals,diff(velocity[2,1,:])./dk,"g-",label="cut g2")
    plot(kvals,diff(diag(velocity[2,:,:]))./dk,"b-",label="cut g2+g2")
    # ylim([-0.02,0.02])
    xlabel("k")
    ylabel("inverse mass")
    legend()
    tight_layout()
    display(fig)
    close(fig)
end
