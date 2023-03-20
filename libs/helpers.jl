using PyPlot 
using CSV
using DataFrames

## run band structure calculations, and save the data 
function get_bandstructure_info(params::Params;ϵ::Float64=0.003,ϕ::Int=0,Da::Float64=4500.0,_hetero::Bool=true,
        fname::String="tmp.csv",lk::Int=110)
    Latt = Lattice()
    initLattice(Latt,params;lk=lk)
    blk = HBM()
    initHBM(blk,Latt,params;lg=9,_σrotation=false)

    df = DataFrame()
    df[!,"kx"] = [real(Latt.kvec[i]) for i in 1:length(Latt.kvec) for j in 1:2]
    df[!,"ky"] = [imag(Latt.kvec[i]) for i in 1:length(Latt.kvec) for j in 1:2]
    df[!,"Hk"] = blk.Hk[:]
    df[!,"Vx"] = blk.Vx[:]
    df[!,"Vy"] = blk.Vy[:]
    CSV.write(fname,df)
    return nothing
end


## plot band structure and velocity related info
function plot_contour_maps(x::Matrix{Float64},y::Matrix{Float64},z::Matrix{Float64};
    contours::Vector{Float64}=[],contour_colors::Vector{String}=["r","g","b"],νs::Vector{Float64}=[1,2,3],
    special_points::Vector{ComplexF64}=[0.0+1im],
    figsize::Tuple=(3,3),_legend::Bool=true,_tight::Bool=true,cmap::String="Spectral",
    savename::String="tmp.pdf")

    fig = figure(figsize=figsize)

    pl = contourf(x,y,z,levels=20,cmap=cmap,vmin=-10,vmax=70)
    if _legend 
        colorbar(pl,ticks=collect(0:10:60))
    end
    
    
    if !isempty(contours)
        pl = contour(x,y,z,levels=contours,colors=contour_colors,linestyles="solid")
        contour_labels = Dict(pl.levels[i]=> @sprintf("%1.1f",νs[i]) for i in eachindex(pl.levels))
        clabel(pl,pl.levels,fmt=contour_labels,inline=true,fontsize=6)
    end
    
    if length(contour_colors)>1
        plot(real(special_points[1:2]),imag(special_points[1:2]),"k*")
        plot(real(special_points[3]),imag(special_points[3]),"o",c=contour_colors[1])
        plot(real(special_points[4]),imag(special_points[4]),"o",c=contour_colors[2])
        plot(real(special_points[5]),imag(special_points[5]),"o",c=contour_colors[3])
        plot(real(special_points[6]),imag(special_points[6]),"kX")
    else
        plot(real(special_points[1:2]),imag(special_points[1:2]),"k*")
        plot(real(special_points[3]),imag(special_points[3]),"o",c="b")
        plot(real(special_points[4]),imag(special_points[4]),"o",c="r")
        plot(real(special_points[5]),imag(special_points[5]),"o",c="g")
        plot(real(special_points[6]),imag(special_points[6]),"kX")
    end

    axis("off")
    axis("equal")
    if _tight
        tight_layout()
    end
    savefig(savename,transparent="true")
    display(fig)
    close(fig)
end


## get special points of a given band structure 
function get_special_points(params::Params;ϵ::Float64=0.003,ϕ::Int=0,Da::Float64=4500.0,_hetero::Bool=true,
    fname::String="tmp.csv",lk::Int=110)
    Latt = Lattice()
    initLattice(Latt,params;lk=lk)
    blk = HBM()
    blk._σrotation=false
    df = DataFrame(CSV.File(fname))
    blk.Hk = reshape(df[!,"Hk"],2,:)
    blk.Vx = reshape(df[!,"Vx"],2,:)
    blk.Vy = reshape(df[!,"Vy"],2,:)
    special_points, special_point_energies = find_special_points_BM(blk,Latt,params)
    special_point_fillings = [[(sum(( sign.(special_point_energies[iband][i] .-blk.Hk) .+1 )./2)/(length(blk.Hk)/2) -1)*4 for i in eachindex(special_point_energies[iband])] for iband in 1:2]
    return special_points, special_point_energies,special_point_fillings
end

## plot special points
function plot_special_points(;fname::String="tmp.csv",flag::String="filling",hline::Float64=20.0, savename::String="tmp.pdf")
    df = DataFrame( CSV.File(fname) )
    fig = figure(figsize=(2.8,2.8))
    colors = ["grey","r","b","g","k","k","g","b","r","grey"]
    iϕs = [5,21]  # strain 003
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
    for i in 1:10
        plot(ϕs,ydata[i,:],"-",c=colors[i],ms=3)
    end 
    # for iϕ in eachindex(ϕs)
    #     if norm(ydata[:,iϕ]) > 1e-5
    #         for i in 1:10
    #             plot(ϕs[iϕ],ydata[i,iϕ],"o",c=colors[i],ms=3)
    #         end
    #     end
    # end
    xlabel(L"$ φ\ (^\circ)$")
    if isequal(flag,"energy")
        ylabel("E (meV)")
    elseif isequal(flag,"filling")
        ylabel("ν")
    end
    if isequal(flag,"energy")
        yticks(-60:30:60)
    else
        ylim([-4.4,4.4])
    end
    # axhline(hline)
    # title("strain004")
    tight_layout()

    savefig(savename,transparent=true)
    display(fig)
    close(fig)

    δ = sort(abs.(ydata[6,:]- ydata[5,:])[:])
    δmax = δ[end]
    δmin = δ[1]
    i = 1
    while δmin < 1e-5
        i = i+1
        δmin = δ[i]
    end
    dνupper = maximum(ydata[9,:]-ydata[8,:])
    dνlower = maximum(ydata[3,:]-ydata[2,:])
    println(δmax,",",δmin,",",dνupper,",",dνlower,",,")
    # println(dνmax)
    return nothing 
end


## plot special points
function plot_special_points_combined(;fname::String="tmp.csv",flag::String="filling",hline::Float64=20.0, savename::String="tmp.pdf")
    df = DataFrame( CSV.File(fname) )
    fig,ax = subplots(2,1,sharex=true,figsize=(2.5,4))
    colors = ["grey","r","b","g","k","k","g","b","r","grey"]
    iϕs = [5,21]  # strain 003
    ϕs = collect(0:2:60)
    ydata = zeros(Float64,10,length(ϕs))
    flags = ["energy","filling"]
    for r in 1:2
        flag = flags[r] 
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
        for i in 1:10
            ax[r].plot(ϕs,ydata[i,:],"-",c=colors[i],ms=3)
        end 
    end

    ax[2].set_xlabel(L"$ φ\ (^\circ)$")
    ax[1].set_ylabel("E (meV)")
    ax[2].set_ylabel("ν")
    ax[1].set_yticks(-60:30:60)
    ax[2].set_ylim([-4.4,4.4])
    tight_layout()

    savefig(savename,transparent=true)
    display(fig)
    close(fig)

    return nothing 
end