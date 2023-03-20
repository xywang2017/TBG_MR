using PyPlot
using Printf
using DataFrames
using CSV
using JLD
using DelimitedFiles
fpath = joinpath(pwd(),"B0")
include(joinpath(fpath,"libs/BM_mod.jl"))
include(joinpath(fpath,"libs/helpers.jl"))

# 
function determine_ν(μ,hk)
    return 4* (sum( (sign.(μ .- hk) .+1 )./2 ) ./(length(hk)/2) -1 )
end
#
lk = 256
ϵs = collect(0.001:0.001:0.007)
ϕs = collect(0:2:60)
Da= -4100.0
iϵ = 2
iϕ = 1
ϵ = ϵs[iϵ]
ϕ = ϕs[iϕ]
str = Int(1000 * ϵ)
str2 = Int(Da)
params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da,_hetero=true)
initParamsWithStrain(params)
Latt = Lattice()
initLattice(Latt,params;lk=lk)
blk = HBM()
fname=joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).csv")
df = DataFrame(CSV.File(fname))
blk.Hk = reshape(df[!,"Hk"],2,:)

fname = joinpath(fpath,"strain/band_structure/eps_$(str)_Da_$(str2)_special_points_$(lk).csv")
df_special_points = DataFrame(CSV.File(fname))
dps = sort(df_special_points[!,"$(ϕ)_k_E"])[6:7]
vhs = sort(df_special_points[!,"$(ϕ)_k_E"])[[4;9;3;10]]
## plot area vs energy curve two bands 
Wmin, Wmax = minimum(blk.Hk), maximum(blk.Hk)
μs = collect(range(0.95Wmin,0.95Wmax,length=150))
νs = zeros(Float64,length(μs))


for i in eachindex(νs)
    νs[i] = determine_ν(μs[i],blk.Hk)
end

fig = figure(figsize=(4,3))
plot(μs,νs,"m-",lw=2)
axvline(dps[1],ls=":",c="k")
axvline(vhs[1],ls=":",c="g")
xlabel("E (meV)")
ylabel("ν")
tight_layout()
display(fig)
close(fig)


###

fig = figure(figsize=(4,3))
plot((μs[1:(end-1)]+μs[2:end])/2,diff(νs)./diff(μs),"m-",lw=2)
axvline(dps[1],ls=":",c="k")
axvline(vhs[1],ls=":",c="g")
xlabel("E (meV)")
ylabel(L"$\mathcal{N}_F$")
tight_layout()
display(fig)
close(fig)


##
νdp_lower = determine_ν(dps[1],blk.Hk)
νvh_lower = determine_ν(vhs[1],blk.Hk)
νdp_upper = determine_ν(dps[2],blk.Hk)
νvh_upper = determine_ν(vhs[2],blk.Hk)
νvh2_lower = determine_ν(vhs[3],blk.Hk)
νvh2_upper = determine_ν(vhs[4],blk.Hk)

## --------------------------------------------------------- ##
fs = load(joinpath(fpath,"strain/transport/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).jld"),"FS")
##
function get_areas(kFS)
    # based on Green's relation ∫ k×dk
    area = 0.0 
    dkFS = diff(kFS)
    kavg = (kFS[1:(end-1)] + kFS[2:end])/2
    for i in eachindex(dkFS)
        area += ( (real(kavg[i])*imag(dkFS[i]) - imag(kavg[i])*real(dkFS[i])) )/2
    end
    area_moire = abs(imag(params.g1'*params.g2))
    return (area) / area_moire  # area in units of moire reciprocal Brillouin zone area
end
function get_FS_areas(fs,νspecial_points)
    # works between -vh1 and vh1 
    areas = []
    for iμ in 1:size(fs,1)
        if fs[iμ][1]>νspecial_points[1] && fs[iμ][1]<νspecial_points[4]
            # fig = figure(figsize=(4,4))
            tmp_area_μ = Float64[]
            for iband in 1:size(fs[iμ][3],1)
                # tmp_area = []
                for iFS in 1:size(fs[iμ][3][iband],1)
                    kFS = fs[iμ][3][iband][iFS][:,1]*params.g1 + fs[iμ][3][iband][iFS][:,2]*params.g2
                    push!(tmp_area_μ,get_areas(kFS)) # based on Green's relation ∫ k×dk
                    # plot(real(kFS)/abs(params.g1),imag(kFS)/abs(params.g1),"b-")
                end
                # push!(tmp_area_μ,copy(tmp_area))
            end
            if length(tmp_area_μ)==2
                push!(areas,[fs[iμ][1], fs[iμ][2],tmp_area_μ])
            else
                push!(areas,[fs[iμ][1], fs[iμ][2],[tmp_area_μ[1],0.0]])
            end
            # tight_layout()
            # display(fig)
            # close(fig)
        end
    end
    return areas
end
areas = get_FS_areas(fs,[νvh_lower,νdp_lower,νdp_upper,νvh_upper])
# for i in 17:25
#     areas[i][3] = areas[i][3][2:-1:1]
# end
for i in eachindex(areas)
    areas[i][3] = sort(areas[i][3])
end

# ----------------------------- 
δν_dirac = 0.0
δν_B = 1/80

ν1_collect = sort([0:δν_B:areas[end][3][1]; -δν_B:-δν_B:areas[1][3][1]])
ν2_collect = sort([0:δν_B:areas[end][3][2]; -δν_B:-δν_B:areas[1][3][2]])
## focus on blue curve, find the corresponding νvalues 
ν2s = [areas[i][3][2] for i in eachindex(areas)]
νs = [areas[i][1] for i in eachindex(areas)]
μs = [areas[i][2] for i in eachindex(areas)]
ν2LLs = []
μ2LLs = []
for i in eachindex(ν2_collect)
    # linear interpolation 
    y = ν2_collect[i]
    j = 1 
    while (y>ν2s[j])
        j = j+1 
    end
    x0,x1, y0,y1 = νs[j-1],νs[j], ν2s[j-1],ν2s[j]
    k = (y1-y0)/(x1-x0)
    b = y0 - k*x0
    push!(ν2LLs,(y-b)/k)

    x0,x1, y0,y1 = μs[j-1],μs[j], ν2s[j-1],ν2s[j]
    k = (y1-y0)/(x1-x0)
    b = y0 - k*x0
    push!(μ2LLs,(y-b)/k)
end

ν1s = [areas[i][3][1] for i in eachindex(areas)]
νs = [areas[i][1] for i in eachindex(areas)]
ν1LLs = []
μ1LLs = []
for i in eachindex(ν1_collect)
    # linear interpolation 
    y = ν1_collect[i]
    j = 1 
    while (y>ν1s[j])
        j = j+1 
    end
    x0,x1, y0,y1 = νs[j-1],νs[j], ν1s[j-1],ν1s[j]
    k = (y1-y0)/(x1-x0)
    b = y0 - k*x0
    push!(ν1LLs,(y-b)/k)

    x0,x1, y0,y1 = μs[j-1],μs[j], ν1s[j-1],ν1s[j]
    k = (y1-y0)/(x1-x0)
    b = y0 - k*x0
    push!(μ1LLs,(y-b)/k)
end

##
fig = figure(figsize=(6,3))
plot([areas[i][1] for i in eachindex(areas)],abs.([areas[i][3][1] for i in eachindex(areas)]),"r.-",label="ν1",ms=4)
plot([areas[i][1] for i in eachindex(areas)],abs.([areas[i][3][2] for i in eachindex(areas)]),"b.-",label="ν2",ms=4)
# plot([areas[i][1] for i in eachindex(areas)],([areas[i][3][1] for i in eachindex(areas)].+[areas[i][3][2] for i in eachindex(areas)]),"g-",label="ν3")
# axvline(νdp_lower,ls=":",c="gray")
# axvline(νdp_upper,ls=":",c="gray")
# axvline(νvh_lower,ls=":",c="gray")
# axvline(νvh_upper,ls=":",c="gray")
axhline(0,ls=":",c="k")
axvline(-0.04)
axvline(0.04)
for i in eachindex(ν1LLs)
    # axvline(ν1LLs[i],ls=":",c="r")
end

for i in eachindex(ν2LLs)
    # axvline(ν2LLs[i],ls=":",c="b")
end
# ylim([-0.02,0.02])
legend()
xlabel("ν")
ylabel("two Fockets")
tight_layout()
savefig("twopockets.pdf",transparent=true)
display(fig)
close(fig)
writedlm("tmp1.txt",[[areas[i][1] for i in eachindex(areas)] [areas[i][3][1]*area_moire/coeff for i in eachindex(areas)] [areas[i][3][2]*area_moire/coeff for i in eachindex(areas)]])
#
# density of states 
ns = [areas[i][1] for i in eachindex(areas)]
μs = [areas[i][2] for i in eachindex(areas)]
n1s = [areas[i][3][1] for i in eachindex(areas)]
n2s = [areas[i][3][2] for i in eachindex(areas)]

ρ1s = diff(n1s) ./ diff(μs)
ρ2s = diff(n2s) ./ diff(μs)
νs_avg = (νs[1:(end-1)] + νs[2:end])/2
fig = figure()
plot(νs_avg,ρ1s,"b-",label=L"ρ_1")
plot(νs_avg,ρ2s,"r-",label=L"ρ_2")
# plot(νs_avg .+ 2νdp_upper,reverse(ρ1s),"r-",label=L"ρ_2")
axvline(0.04)
axvline(-0.04)
legend()
xlabel("ν")
ylabel("DoS")
tight_layout()
display(fig)
close(fig)

## -------------------- 
function get_FS_areas_pastvH1(fs,νspecial_points)
    # works between -vh1 and vh1 
    areas = []
    for iμ in 1:size(fs,1)
        if fs[iμ][1]<νspecial_points[1] && fs[iμ][1]>νspecial_points[5]
            # fig = figure(figsize=(4,4))
            tmp_area_μ = Float64[]
            for iband in 1:size(fs[iμ][3],1)
                # tmp_area = []
                for iFS in 1:size(fs[iμ][3][iband],1)
                    kFS = fs[iμ][3][iband][iFS][:,1]*params.g1 + fs[iμ][3][iband][iFS][:,2]*params.g2
                    push!(tmp_area_μ,get_areas(kFS)) # based on Green's relation ∫ k×dk
                    # plot(real(kFS)/abs(params.g1),imag(kFS)/abs(params.g1),"b-")
                end
                # push!(tmp_area_μ,copy(tmp_area))
            end
            if length(tmp_area_μ)==1
                push!(areas,[fs[iμ][1], fs[iμ][2],tmp_area_μ])
            else
                println("error with FS count")
            end
            # tight_layout()
            # display(fig)
            # close(fig)
        end
    end
    return areas
end
areas = get_FS_areas_pastvH1(fs,[νvh_lower,νdp_lower,νdp_upper,νvh_upper,νvh2_lower,νvh2_upper])
ħ= 1.054571817e-34 
q = 1.6e-19
area_moire = imag(params.g1'*params.g2)/(2.46e-10)^2
coeff = 2π*q/ħ 
writedlm("tmp3.txt",[[areas[i][1] for i in eachindex(areas)] [areas[i][3][1]*area_moire/coeff for i in eachindex(areas)] ])
##

fig = figure()
data = readdlm("tmp2.txt")
plot(data[:,1],data[:,2])
data = readdlm("tmp3.txt")
plot(data[:,1],abs.(data[:,2]))
data = readdlm("tmp1.txt")
plot(data[:,1],abs.(data[:,2]))
plot(data[:,1],abs.(data[:,3]))
plot(data[:,1],abs.(data[:,3].+data[:,2]))
display(fig)
close(fig)
ν3_collect = collect((areas[1][3][1]+1e-3):δν_B:areas[end][3][1])
## focus on blue curve, find the corresponding νvalues 
ν3s = [areas[i][3][1] for i in eachindex(areas)]
νs = [areas[i][1] for i in eachindex(areas)]
μs = [areas[i][2] for i in eachindex(areas)]
ν3LLs = []
μ3LLs = []
for i in eachindex(ν3_collect)
    # linear interpolation 
    y = ν3_collect[i]
    j = 1 
    while (y>ν3s[j])
        j = j+1 
    end
    x0,x1, y0,y1 = νs[j-1],νs[j], ν3s[j-1],ν3s[j]
    k = (y1-y0)/(x1-x0)
    b = y0 - k*x0
    push!(ν3LLs,(y-b)/k)

    x0,x1, y0,y1 = μs[j-1],μs[j], ν3s[j-1],ν3s[j]
    k = (y1-y0)/(x1-x0)
    b = y0 - k*x0
    push!(μ3LLs,(y-b)/k)
end

##
sort_indx = sortperm(Float64[ν1LLs;ν2LLs;ν3LLs])
combined_νLLs = sort(Float64[ν1LLs;ν2LLs;ν3LLs])
combined_μLLs = Float64[μ1LLs;μ2LLs;μ3LLs][sort_indx]



### -----------------------------------------
# xE = range(0.9*vhs[1],0.99*vhs[4],length=600)
# xν = [determine_ν(xE[i],blk.Hk) for i in eachindex(xE)]
# δ = 1
# y = [ sum(1/π * δ./((xE[i].-combined_μLLs).^2 .+δ^2)) for i in eachindex(xE)]
# # y = [ sum(1/π * δ./((xν[i].-combined_νLLs).^2 .+δ^2)) for i in eachindex(xE)]

# fig,ax = subplots(2,1,sharex=true,figsize=((8,6)))
# ax[1].plot(xν, 1 ./ y,"b-")
# ax[1].set_yscale("log")
# max_i = floor(νvh2_upper/(4δν_B))
# for i in -max_i:max_i 
#     ax[1].axvline(i*(4δν_B),c="r",ls=":")
# end
# ax[1].axvline(νvh_upper,c="g",ls="-")
# ax[1].axvline(νdp_upper,c="g",ls="-")
# ax[1].axvline(νdp_lower,c="g",ls="-")
# # ax[1].set_xlabel("ν")
# ax[1].set_ylabel("")
# ax[1].set_ylim([0,0.04])


# for i in eachindex(combined_μLLs)
#     # ax[2].axhline(combined_μLLs[i],ls="-",c="gray")
#     # ax[2].axvline(combined_νLLs[i],ls="-",c="gray")
#     ax[2].plot(combined_νLLs[i],combined_μLLs[i],"b+",ms=6,lw=3)
# end

# ax[2].plot(νdp_lower,dps[1],"b+",lw=3)
# ax[2].axhline(dps[2],ls="-",c="b")
# ax[2].plot(xν,xE,"k-")
# ax[2].set_xlabel("ν")
# ax[2].set_ylabel("E (meV)")
# ax[2].set_xlim([-0.6,1.4])
# ax[2].set_ylim([-11,18])
# tight_layout()
# display(fig)
# close(fig)


fig = figure(figsize=(2,7))
data = readdlm("data_Strain/p1q80_strain_003.txt")
energies = reshape(data[:,1],160,:)
plot(ones(length(energies)),energies[:],"r_")
# plot(ones(length(combined_μLLs)),combined_μLLs,"b_",ms=6,lw=6)
axhline(vhs[1])
axhline(vhs[2])
ylim([-20,20])
ylabel("E (meV)")
display(fig)
savefig("tmp1.pdf",transparet=true)
close(fig)

fig = figure(figsize=(6,3))
data = readdlm("data_Strain/p1q80_strain_003.txt")
energies = reshape(data[:,1],160,:)
axhline(vhs[1])
axhline(vhs[2])
plot(eachindex(energies[:,1]),energies[:,1],"r.",ms=1)
# plot(ones(length(combined_μLLs)),combined_μLLs,"b_",ms=6,lw=6)

# ylim([-20,20])
# ylabel("E (meV)")
display(fig)
savefig("tmp.pdf",transparet=true)
close(fig)