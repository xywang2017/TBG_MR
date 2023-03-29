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
iϵ = 2
ϵ = ϵs[iϵ]
str = Int(1000 * ϵ)
str2 = Int(Da)

function get_dos(energy,hk,δ)
    dos = 1/π * sum( δ ./ ((energy .- hk).^2 .+ δ^2) )/ (length(hk)/2)
    return dos
end

function get_filling(energy,hk)
    ν = (sum((sign.(energy .- hk) .+ 1)./2) / (length(hk)/2) -1)*4
end

function dos_map(ϕs,δ)
    energies = range(-70,70,length=400)
    νs = zeros(Float64,length(energies),length(ϕs))
    dos = zeros(Float64,length(energies),length(ϕs))
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        fname = joinpath(fpath,"strain/band_structure/w0_0_eps_$(str)_phi_$(ϕ)_Da_$(str2).csv")
        df = DataFrame(CSV.File(fname))
        hk = df[!,"Hk"]
        for j in eachindex(energies)
            dos[j,iϕ] = get_dos(energies[j],hk,δ)
            νs[j,iϕ] = get_filling(energies[j],hk)
        end
    end

    fname = joinpath(fpath,"strain/band_structure/w0_0_eps_$(str)_Da_-4100_special_points.csv")
    df = DataFrame( CSV.File(fname) )
    flag = "energy"
    iϕs = [6,22]  # strain 002
    iϕs = [5,21]  # strain 003
    iϕs = [6,22]  # strain 003
    # iϕs = [20,28]  # strain 007
    ϕs = collect(0:2:60)
    ydata = zeros(Float64,10,length(ϕs))
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        if "$(ϕ)_k_E" in names(df)
            if isequal(flag,"energy")
                data = sort(df[!,"$(ϕ)_k_E"])
            elseif isequal(flag,"filling")
                data = sort(df[!,"$(ϕ)_k_nu"])
            end
            y = []
            if iϕ <= iϕs[1] 
                push!(y,data[1])
                push!(y,data[2])
                push!(y,data[3])
                push!(y,data[4])
                push!(y,(data[5]+data[6])/2)
                push!(y,(data[7]+data[8])/2)
                push!(y,data[9])
                push!(y,data[10])
                push!(y,data[11])
                push!(y,data[12])
            elseif (iϕ > iϕs[1] && iϕ<=iϕs[2])
                push!(y,data[1])
                push!(y,data[2])
                push!(y,data[4])
                push!(y,data[3])
                push!(y,(data[5]+data[6])/2)
                push!(y,(data[7]+data[8])/2)
                push!(y,data[10])
                push!(y,data[9])
                push!(y,data[11])
                push!(y,data[12])
            elseif (iϕ>=iϕs[2]) 
                push!(y,data[1])
                push!(y,data[3])
                push!(y,data[4])
                push!(y,data[2])
                push!(y,(data[5]+data[6])/2)
                push!(y,(data[7]+data[8])/2)
                push!(y,data[11])
                push!(y,data[9])
                push!(y,data[10])
                push!(y,data[12])
            end
            ydata[:,iϕ] = y
        end
    end
    fig = figure(figsize=(2.8,2.8))
    pcolormesh(ϕs,energies,dos,linewidth=0,rasterized=true)
    # pcolormesh(repeat(ϕs',length(energies)),νs,dos,linewidth=0,rasterized=true)
    colors = ["white","r","b","g","k","k","g","b","r","white"]
    for i in 1:10
        # plot(ϕs[ydata[i,:].!=0],ydata[i,ydata[i,:].!=0],".",c=colors[i],ms=3,lw=1)
        plot(ϕs,ydata[i,:],"--",c=colors[i],ms=3,lw=1)
    end 
    # axhline(1.4)
    # axhline(0.6)
    # axhline(0.1)
    # axvline(0)
    xlabel(L"$\varphi\ (^\circ)$")
    ylabel("ν")
    # ylim([0,3])
    # axhline(0.8,ls=":",c="k")
    # axhline(1.7,ls=":",c="k")
    # axhline(2.5,ls=":",c="k")
    # yticks([0.8,1.7,2.5])
    # colorbar()
    # ylabel("ν")
    # ylim([-2,2])
    title("ϵ=0.2%")
    tight_layout()
    # savefig("dos_color_plot_filling_004.pdf",transparent=true)
    # savefig("dos_color_plot_filling_004.png",transparent=true,dpi=600)
    display(fig)
    close(fig)
    
end

dos_map(ϕs,0.3)


function dos_mapv1(ϕs,δ)
    energies = range(-70,70,length=400)
    νs = zeros(Float64,length(energies),length(ϕs))
    dos = zeros(Float64,length(energies),length(ϕs))
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        fname = joinpath(fpath,"strain/band_structure/w0_0_eps_$(str)_phi_$(ϕ)_Da_$(str2).csv")
        df = DataFrame(CSV.File(fname))
        hk = df[!,"Hk"]
        for j in eachindex(energies)
            dos[j,iϕ] = get_dos(energies[j],hk,δ)
            νs[j,iϕ] = get_filling(energies[j],hk)
        end
    end

    fig = figure(figsize=(2.8,2.8))
    # pcolormesh(ϕs,energies,dos,linewidth=0,rasterized=true)
    pcolormesh(repeat(ϕs',length(energies)),νs,dos,linewidth=0,rasterized=true)
    xlabel(L"$\varphi\ (^\circ)$")
    ylabel("ν")
    ylim([0,3])
    axhline(0.8,ls=":",c="k")
    axhline(1.7,ls=":",c="k")
    axhline(2.5,ls=":",c="k")
    yticks([0.8,1.7,2.5])
    title("ϵ=0.2%")
    tight_layout()
    # savefig("dos_color_plot_filling_004.pdf",transparent=true)
    savefig("dos_color_plot_filling_002_w00.png",transparent=true,dpi=600)
    display(fig)
    close(fig)
    
end

dos_mapv1(ϕs,0.3)