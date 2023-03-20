using PyPlot
using Printf
using DataFrames
using CSV
using DelimitedFiles
using JLD
# using Peaks
fpath = joinpath(pwd(),"B0")
include(joinpath(fpath,"libs/BM_mod.jl"))
include(joinpath(fpath,"libs/helpers_transport.jl"))
include(joinpath(fpath,"libs/helpers.jl"))
# basic settings
ϵ, ϕ, Da = 0.002, 0, -4100.0
str, str2 = Int(1000*ϵ), Int(Da)
lk = 256
params = Params(ϵ=ϵ,φ=ϕ*π/180,Da=Da,dθ=1.38π/180)
initParamsWithStrain(params)
Latt = Lattice()
initLattice(Latt,params;lk=lk)
blk = HBM()
# initHBM(blk,Latt,params;lg=9,_σrotation=false)
# fname=joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).csv")
fname=joinpath(fpath,"strain/band_structure/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).csv")

# df = DataFrame()
# df[!,"kx"] = [real(Latt.kvec[i]) for i in 1:length(Latt.kvec) for j in 1:2]
# df[!,"ky"] = [imag(Latt.kvec[i]) for i in 1:length(Latt.kvec) for j in 1:2]
# df[!,"Hk"] = blk.Hk[:]
# df[!,"Vx"] = blk.Vx[:]
# df[!,"Vy"] = blk.Vy[:]
# CSV.write(fname,df)
df = DataFrame(CSV.File(fname))
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
eBτs = [collect(-10:1:-2);collect(-1:0.1:-0.1);collect(0.1:0.1:1);collect(2:1:10)].*1e-4
transport_coefficients,FS_collect = mainTransport(blk,Latt,params;eBτs = eBτs,nμs=121);
save(joinpath(fpath,"strain/transport/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).jld"),"transport",transport_coefficients,"FS",FS_collect)
# save(joinpath(fpath,"strain/transport/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).jld"),"transport",transport_coefficients)
ϵa = 1260.88
ω0τs = eBτs * ϵa
transport_coefficients = load(joinpath(fpath,"strain/transport/eps_$(str)_phi_$(ϕ)_Da_$(str2)_$(lk).jld"),"transport")

### Plots #### 
fname_special_points = joinpath(fpath,"strain/band_structure/eps_$(str)_Da_$(str2)_special_points.csv")

## plot velocity fields 
fig = figure(figsize=(4,3))
colors=["tab:red","tab:blue","tab:green"]
cnt = 1
data1 = load("closedFS3.48.jld")
vx = data1["vx"][1]
vy = data1["vy"][1]
s = data1["s"][1]
plot(s./s[end],vx/1000,c=colors[cnt],label=L"\rm closed\ FS")
plot((s[1:(end-1)]+s[2:end])./(2s[end]),ones(length(vx)-1)*sum(0.5*(vx[1:(end-1)]+vx[2:end]).*diff(s))/(s[end])/1000,":",c=colors[cnt])
cnt+=1
data2 = load("openFS2.47.jld")
for i in 1:size(data2["vx"],1)
        vx = data2["vx"][i]
        vy = data2["vy"][i]
        s = data2["s"][i]
        if i==1 
                label = L"\rm open\ FS1"
        else 
                label = L"\rm open\ FS2"
        end
        plot(s./s[end],vx/1000,c=colors[cnt],label=label)
        plot((s[1:(end-1)]+s[2:end])./(2s[end]),ones(length(vx)-1)*sum(0.5*(vx[1:(end-1)]+vx[2:end]).*diff(s))/(s[end])/1000,":",c=colors[cnt])
        cnt+=1
end
axvline(1,ls="-",c="k",lw=1)
xlabel(L"s/s_0")
xticks([0,1],["0",L"1"])
xlim([0,1.75])
ylabel(L"\rm ħv_x/a\ (eV)")
legend(loc="upper right")
tight_layout()
savefig("velocity_open_vs_closedFS.pdf",transparent=true)
display(fig)
close(fig)

###
fig = figure(figsize=(4,4))
kvec = real(Latt.kvec)*params.g1 .+ imag(Latt.kvec)*params.g2
x = real(reshape(kvec,Latt.lk,Latt.lk))
y = imag(reshape(kvec,Latt.lk,Latt.lk))
z = reshape(blk.Hk[1,:],Latt.lk,Latt.lk)
colors=["tab:red","tab:blue","tab:green"]
cnt = 1
for i in 1:size(data1["FS"],1)
        plot(real(data1["FS"][i]),imag(data1["FS"][i]),".",c=colors[cnt],ms=2)
        cnt+=1
end
for i in 1:size(data2["FS"],1)
        plot(real(data2["FS"][i]),imag(data2["FS"][i]),".",c=colors[cnt],ms=2)
        cnt +=1
end

plot([0,real(params.g1)],[0,imag(params.g1)],"k-")
plot([0,real(params.g2)],[0,imag(params.g2)],"k-")
plot(real(params.g1).+[0,real(params.g2)],imag(params.g1).+[0,imag(params.g2)],"k-")
plot(real(params.g2).+[0,real(params.g1)],imag(params.g2).+[0,imag(params.g1)],"k-")
axis("equal")
axis("off")
tight_layout()
savefig("contours.pdf",transparent=true)
display(fig)
close(fig)

# unsymmetrized transport
plot_ρs(transport_coefficients,eBτs,fname_special_points)

# symmetrization
νs,ρ,δρ,δθ = Bsymmetrization(transport_coefficients,eBτs) ;

# symmetrized transport
plot_symmetrized_transport(νs,ρ,fname_special_points)

# hall
plot_hall(νs,δρ,fname_special_points)

# ωcτ
plot_ωcτavg_versus_filling(νs,fname_special_points)

# principal axis
plot_principle_strain_axis(νs,δθ,fname_special_points)

# geometric --- principal axis plotted against moire unit cell vectors
plot_principal_strain_axis_vs_unit_cell(110)

# line cuts at a fixed filling, plot (symmetrized/unsymmetrized) rho vs B
plot_transport_filling(νs,ρ,"quadratic",2)

# multiple line cuts
# plot_transport_filling(νs,ρ,"quadratic",52:2:64)