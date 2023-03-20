# plot ρxx, ρyy vs field eBτ
function plot_ρs(transport_coefficients,eBτs,fname)
    fig,ax = subplots(figsize=(5,3))
    pl1, pl2 = 0,0
    colors1 = [[0.;0.6;0.]*iτ/(length(eBτs)/2) for iτ in 1:((length(eBτs)÷2))]
    colors2 = [[0.8;0;0.8]*iτ/(length(eBτs)/2) for iτ in 1:((length(eBτs)÷2))]
    colors1 = reverse(["tab:blue","tab:red","tab:green","tab:purple"])
    # colors1 = reverse(["blue","red","green","purple"])
    νs = [4transport_coefficients[i][1]-4 for i in 1:size(transport_coefficients,1)]
    
    idx_select = [6;9;11;18]
    cnt = 1
    for iτ in reverse(idx_select)
        str = @sprintf("%.2f",abs(ω0τs[iτ]))
        ρxx = [transport_coefficients[i][2][1,1,iτ] for i in 1:size(transport_coefficients,1)]
        ρyy = [transport_coefficients[i][2][2,2,iτ] for i in 1:size(transport_coefficients,1)]
        pl1, = ax.plot(νs,ρxx,":",c=colors1[cnt],ms=1,lw=2)
        pl2, = ax.plot(νs,ρyy,"-",c=colors1[cnt],label=str,ms=1,lw=2)
        cnt = cnt + 1
    end
    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    # data = sort(df[!,"$(ϕ)_k_nu"])[[2;9]]
    for i in 1:length(data)
        ax.axvline(data[i],ls=":",c="gray")
    end
    ax.legend(loc="upper left",title=L"ω_0τ=")
    # ax.legend(loc="upper left",title=L"δB/B=")
    ax.set_xlabel("ν")
    ax.set_ylabel(L"$ρ/ \left(ρ_Q \frac{Γ}{ϵ_M} \right)$",rotation=90)
    ax.set_yscale("log")
    ax.set_ylim([1.2,0.8e4])
    # ax.set_ylim([0,80])
    # ax.set_ylim([2,0.8e2])
    # subplots_adjust(wspace=0, hspace=0)
    tight_layout()
    # savefig("strained_unsymmetrized_rho_phi4.pdf",transparent=true)
    # savefig("unstrained_unsymmetrized_rho.pdf",transparent=true)
    display(fig)
    close(fig)
    return nothing
end


## symmetrization wrt B to -B
function Bsymmetrization(transport_coefficients,eBτs)
    ρ = zeros(2,2,length(eBτs)÷2,size(transport_coefficients,1))
    δρ = zeros(2,2,length(eBτs)÷2,size(transport_coefficients,1))
    νs = zeros(size(transport_coefficients,1))
    for iμ in 1:size(transport_coefficients,1)
        νs[iμ] = transport_coefficients[iμ][1]
        for iτ in 1:(length(eBτs)÷2)
            ρ[:,:,iτ,iμ] = (transport_coefficients[iμ][2][:,:,iτ] + transport_coefficients[iμ][2][:,:,length(eBτs)-iτ+1] )/2
            δρ[:,:,iτ,iμ] = (transport_coefficients[iμ][2][:,:,length(eBτs)-iτ+1] -transport_coefficients[iμ][2][:,:,iτ] )/2 
        end
    end
    # principle axis
    Uθ = zeros(2,2,length(eBτs)÷2,size(transport_coefficients,1))
    δθ = zeros(length(eBτs)÷2,size(transport_coefficients,1))
    for iμ in 1:size(transport_coefficients,1)
        for iτ in 1:(length(eBτs)÷2)
            if (abs(ρ[1,2,iτ,iμ]-ρ[2,1,iτ,iμ])>1e-6)
                println("error with symmetrization ρ[1,2]≂̸ρ[2,1]")
            end
            if (abs(δρ[1,1,iτ,iμ])+abs(δρ[2,2,iτ,iμ])>1e-6)
                println("error with antisymmetrization δρ[1,1]≂̸δρ[2,2]≂̸0")
            end
            if (abs(δρ[1,2,iτ,iμ]+δρ[2,1,iτ,iμ])>1e-6)
                println("error with antisymmetrization δρ[1,2]≂̸-δρ[2,1]")
            end
            ρ11,ρ12,ρ22 = ρ[1,1,iτ,iμ],ρ[1,2,iτ,iμ],ρ[2,2,iτ,iμ]
            ω = sqrt((ρ11-ρ22)^2/4+ρ12^2)
            ϵminus = (ρ11+ρ22)/2 - ω
            # δθ[iτ,iμ] = - atan((ρ11-ϵminus)/ρ12)  # (-π/2,π/2)
            δθ[iτ,iμ] = acos(-sign(ρ12)*sqrt(0.5*(1-(ρ11-ρ22)/(2ω))))
            Uθ[:,:,iτ,iμ] = Float64[cos(δθ[iτ,iμ]) -sin(δθ[iτ,iμ]);sin(δθ[iτ,iμ]) cos(δθ[iτ,iμ])]
            ρ[:,:,iτ,iμ] = Uθ[:,:,iτ,iμ]' * ρ[:,:,iτ,iμ] * Uθ[:,:,iτ,iμ]
            # δρ[:,:,iτ,iμ] = Uθ[:,:,iτ,iμ]' * δρ[:,:,iτ,iμ] * Uθ[:,:,iτ,iμ]
            
            if (abs(ρ[1,2,iτ,iμ])+abs(ρ[2,1,iτ,iμ])) > 1e-5
                println("error with diagonalization using U")
            end
            if (abs(δρ[1,1,iτ,iμ])+abs(δρ[2,2,iτ,iμ])) > 1e-5
                println("error with off-diagonalization using U")
            end
            if abs(δρ[1,2,iτ,iμ]+δρ[2,1,iτ,iμ]) > 1e-5 
                println("error with antisymmetrization using U")
            end
        end
    end
    return νs,ρ,δρ,δθ
end


function plot_symmetrized_transport(νs,ρ,fname)
    fig,ax = subplots(figsize=(5,3))
    pl1, pl2 = 0,0
    colors1 = [[1.;0.0;0.]*iτ/(length(eBτs)/2) for iτ in 1:((length(eBτs)÷2))]
    colors2 = [[0;0.0;1]*iτ/(length(eBτs)/2) for iτ in 1:((length(eBτs)÷2))]
    colors1 = reverse(["tab:blue","tab:red","tab:green","tab:purple"])
    # for iτ in 1:2:((length(eBτs)÷2))
    idx_select = [6;9;11;18]
    # idx_select = [6;9;11;18]
    cnt = 1
    for iτ in reverse(idx_select)
        str = @sprintf("%.2f",-ω0τs[iτ])
        pl1, = ax.plot(4νs.-4,ρ[1,1,iτ,:],":",c=colors1[cnt],ms=3,lw=2)
        pl2, = ax.plot(4νs.-4,ρ[2,2,iτ,:],"-",c=colors1[cnt],ms=3,label=str,lw=2)
        cnt += 1
    end

    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    for i in 1:length(data)
       ax.axvline(data[i],ls=":",c="gray")
    end

    ax.set_xlabel("ν")
    ax.set_ylabel(L"$ρ/ \left(ρ_Q \frac{Γ}{ϵ_M} \right)$")
    
    ax.legend(loc="upper left",title=L"ω_0τ=")
    ax.set_ylim([1.2,0.8e4])
    ax.set_yscale("log")
    # xlim([-1,1])
    tight_layout()
    # savefig("symmetrize_rho.pdf",transparent=true)
    display(fig)
    close(fig)
end


##
function plot_hall(νs,δρ,fname)
    fig = figure(figsize=(5,3))
    colors = [[(size(ρ,3)-iτ)/(size(ρ,3)-1);0.3;(iτ-1)/(size(ρ,3)-1)] for iτ in 1:(size(ρ,3))]
    colors = reverse(["tab:blue","tab:red","tab:green","tab:purple"])
    # colors = [[0.8;0.2;0.2]*(iτ-1)/(size(ρ,3)-1) for iτ in 1:(size(ρ,3))]
    pl1, pl2 = 0,0
    area = abs(imag(params.a1'*params.a2)) 
    plot(4νs.-4,(4νs.-4),"--",c="k")
    plot(4νs.-4,4νs,"--",c="k")
    plot(4νs.-4,4νs.-8,"--",c="k")
    # for iτ in 1:3:size(ρ,3)
    idx_select = [6;9;11;18]
    # idx_select = [7;8;9;10]
    cnt = 1
    for iτ in reverse(idx_select)
    # for iτ in [1;4;8;9;10;collect(12:3:17)]
        str = @sprintf("%.2f",-ω0τs[iτ])
        nH = - eBτs[iτ] ./ ( δρ[2,1,iτ,:] .*( 2π / (params.vf*params.kb) ) ) * area # /(area * 2.46^2*1e-16)
        # nH =  -δρ[2,1,iτ,:]
        # idx_toplot = [collect(1:(length(nH)÷2));collect((length(nH)÷2+2):length(nH))]
        # pl1, = plot(4νs[idx_toplot].-4, nH[idx_toplot] ,"-",c=colors[iτ],ms=3,label=str)
        pl1, = plot(4νs.-4, nH ,"-",c=colors[cnt],ms=1,label=str,lw=2)
        cnt += 1
    end
    
    # axhline(0,c="k",ls="-",lw=0.5)
    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    # data = sort(df[!,"$(ϕ)_k_nu"])[[2;9]]
    for i in 1:length(data)
        axvline(data[i],ls=":",c="gray")
    end
    legend(loc="upper left",title=L"ω_0τ=")
    xlabel("ν")
    ylabel(L"$ρ_H/ \left(ρ_Q \frac{Γ}{ϵ_M} \right)$")
    # ylabel(L"\rm n_H\ (moir\'e\ u.c.)")
    # ylabel(L" ρ_{21}")
    ylim([-6,6])
    # xlim([-4,4])
    tight_layout()
    # savefig("strained_rho_21.pdf",transparent=true)
    display(fig)
    close(fig)
end


##
function plot_ωcτavg_versus_filling(νs,fname)
    fig = figure(figsize=(5,3))
    colors = [[cos(π/2* iτ/(size(ρ,3)));0.4;sin(π/2* iτ/(size(ρ,3)))] for iτ in 1:(size(ρ,3))]
    ωcτs_avg = zeros(length(νs),length(eBτs))
    for id in eachindex(νs) 
        # println(id)
        # if ! (id in [60;61])
        ωcτ = reshape(transport_coefficients[id][3][1],:,length(eBτs))
        nFS = size(ωcτ,1)
        ωcτs_avg[id,:] = reshape(sum(abs.(ωcτ),dims=1) / nFS ,:)
        # end
    end
    idx_select = [6]
    for iτ in idx_select
        str = @sprintf("ω0τ=%.2f",-ω0τs[iτ])
        plot(4νs.-4,-ωcτs_avg[:,iτ]./ω0τs[iτ],"-",c="k",ms=3,lw=2)
    end
    # yticks([])

    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    # data = sort(df[!,"$(ϕ)_k_nu"])[[2;9]]
    for i in 1:length(data)
        axvline(data[i],ls=":",c="gray")
    end
    ylim([0,20])
    # legend()
    ylabel(L"m_0/m^*")
    xlabel(L"ν")
    tight_layout()
    savefig("omegactau.pdf",transparent=true)
    display(fig)
    close(fig)
end


##
function plot_principle_strain_axis(νs,δθ,fname)
    fig = figure(figsize=(5,3))
    l0,l1,l2,l3=0,0,0,0
    colors = [[(size(ρ,3)-iτ)/(size(ρ,3));0.3;iτ/(size(ρ,3))] for iτ in 1:(size(ρ,3))]
    colors = reverse(["tab:blue","tab:red","tab:green","tab:purple"])
    idx_select = [6;9;11;19]
    cnt = 1
    for iτ in reverse(idx_select) #1:(size(ρ,3))
    # iτ = 3
        str = @sprintf("%.2f",-ω0τs[iτ])
        # if iτ != 6
            plot(4νs[3:(end-3)].-4,δθ[iτ,3:(end-3)]*180/π,"-",c=colors[cnt],ms=5,label=str,lw=2)
        # else
        #     plot(4νs[abs.(4νs.-4).>0.06].-4,δθ[iτ,abs.(4νs.-4).>0.06]*180/π,"-",c=colors[cnt],ms=5,label=str,lw=2)
        # end
        cnt += 1
    end
    
    l1,=plot(4νs.-4,ones(length(νs))*angle(params.a1)*180/π,"--")
    # l1,=plot(4νs.-4,ones(length(νs))*(90+angle(params.a1)*180/π),"--")

    l2,=plot(4νs.-4,ones(length(νs))*angle(params.a2)*180/π,"--")
    l3,=plot(4νs.-4,ones(length(νs))*angle(params.a2-params.a1)*180/π,"--")
    # plot(4νs.-4,ones(length(νs))*180,"--")
    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    for i in 1:length(data)
        axvline(data[i],ls=":",c="gray")
    end
    legend(loc="upper right",title=L"ω_0τ=")
    # legend([l1,l2,l3],[L"a_1",L"a_2",L"a_3"],loc="upper left")
    xlabel("ν")
    ylabel("δθ")
    yticks([0,60,120],[L"$0$",L"$\frac{π}{3}$",L"$\frac{2π}{3}$"])
    ylim([0,185])
    tight_layout()
    savefig("principal_axis.pdf",transparent=true)
    display(fig)
    close(fig)
end


function plot_principal_strain_axis_vs_unit_cell(iν::Int)
    fig = figure(figsize=(2.5,2.5))
    dθ = 1.38π/180
    kb = 8π/3*sin(dθ/2)
    a1 = 4π/(3kb)*exp(1im * π/6)
    a2 = 4π/(3kb)*1im

    lattice_points0 = [0.0+0.0im,a1,a2,a2-a1,-a1,-a2,a1-a2]
    lattice_points = [0.0+0.0im,params.a1,params.a2,params.a2-params.a1,
                        -params.a1,-params.a2,params.a1-params.a2]

    scatter(real(lattice_points0),imag(lattice_points0),s=20,c="grey")
    plot(real([lattice_points0[1],lattice_points0[2]]),imag([lattice_points0[1],lattice_points0[2]]),"--",c="grey")
    plot(real([lattice_points0[1],lattice_points0[3]]),imag([lattice_points0[1],lattice_points0[3]]),"--",c="grey")
    plot(real([lattice_points0[3],lattice_points0[2]]),imag([lattice_points0[3],lattice_points0[2]]),"--",c="grey")

    scatter(real(lattice_points),imag(lattice_points),s=20,c="k")
    plot(real([lattice_points[1],lattice_points[2]]),imag([lattice_points[1],lattice_points[2]]),"k-")
    plot(real([lattice_points[1],lattice_points[3]]),imag([lattice_points[1],lattice_points[3]]),"k-")
    plot(real([lattice_points[3],lattice_points[2]]),imag([lattice_points[3],lattice_points[2]]),"k-")

    # plot(0.5*real([lattice_points0[5],lattice_points0[7]]),imag([lattice_points0[5],lattice_points0[7]]),"-",c="k")
    # text(0,-18,"δθ")
    ν = 4νs[iν]-4
    θ0 = δθ[end,iν]
    radius = abs(params.a1)*0.7
    plot(radius*[-cos(θ0),cos(θ0)],radius*[-sin(θ0),sin(θ0)],"r-")
    plot(radius*[-sin(θ0),sin(θ0)],radius*[cos(θ0),-cos(θ0)],"b-")

    axis("off")
    str = @sprintf "%.2f" ν
    title(L"$ϕ=%$(ϕ)^\circ,\nu\approx %$(str)$")
    title(L"$\nu\approx %$(str)$")
    axis("equal")
    tight_layout()
    savefig("eps02_$(ν).pdf",transparent=true)
    display(fig)
    close(fig)
    return nothing 
end

## versus field plot at a particular filling
function plot_transport_filling(νs,ρ,flag,ν)
    fig = figure(figsize=(3,3))
    id = argmin(abs.(4νs .-4 .- ν)) 
    if size(transport_coefficients[id][3],1) == 2 
        println("semimetal")
    end
    ωcτs = reshape(transport_coefficients[id][3][1],:,length(eBτs))
    nFS = size(ωcτs,1)
    ωcτ_avg = reshape(sum(ωcτs,dims=1) / nFS ,:)
    
    ρ1 = ρ[1,1,:,id]
    ρ2 = ρ[2,2,:,id]
    ρ1 = transport_coefficients[id][2][1,1,1:length(eBτs)÷2]
    ρ2 = transport_coefficients[id][2][2,2,1:length(eBτs)÷2]


    # versus eBτ a^2/ħ^2
    if isequal(flag,"versusB")
        pl1, = plot(abs2.(eBτs[1:length(eBτs)÷2])*1e8,ρ1,"r-^",ms=4)
        pl2, = plot(abs2.(eBτs[1:length(eBτs)÷2])*1e8,ρ2,"b->",ms=4)
        xlabel(L"$\left(eBτ\frac{a^2}{ħ^2}\times 10^4\right)^2$")
    elseif isequal(flag,"quadratic")
        pl1, = plot(abs2.(ωcτ_avg[1:length(eBτs)÷2]),ρ1,"r-",ms=4,lw=2)
        pl2, = plot(abs2.(ωcτ_avg[1:length(eBτs)÷2]),ρ2,"b-",ms=4,lw=2)
        xlabel(L"(ω̄_cτ)^2")
    elseif isequal(flag,"linear")
        pl1, = plot(abs.(ωcτ_avg[1:length(eBτs)÷2]),ρ1,"r-^",ms=4)
        pl2, = plot(abs.(ωcτ_avg[1:length(eBτs)÷2]),ρ2,"b->",ms=4)
        xlabel(L"ω̄_cτ")
    end
    
    # yscale("log")
    # xscale("log")
    ylabel(L"$ρ/ \left(ρ_Q \frac{Γ}{ϵ_M} \right)$")
    # legend([pl1,pl2],[L"$ρ_{1}$",L"$ρ_{2}$"])
    legend([pl1,pl2],[L"$ρ_{xx}$",L"$ρ_{yy}$"])
    title("ν=$(ν)")
    tight_layout()
    savefig("nu0.5_unsymmetrized.pdf",transparent=true)
    display(fig)
    close(fig)
end


## versus field plot for a range of filling
function plot_transport_filling_range(νs,ρ,flag,iνs)
    fig = figure(figsize=(4,3))
    for iν in iνs 
        id = argmin(abs.(νs .- νs[iν])) 
        if size(transport_coefficients[id][3],1) == 2 
            println("semimetal")
        end
        ωcτs = reshape(transport_coefficients[id][3][1],:,length(eBτs))
        nFS = size(ωcτs,1)
        ωcτ_avg = reshape(sum(ωcτs,dims=1) / nFS ,:)
        
        ρ1 = ρ[1,1,:,id]
        ρ2 = ρ[2,2,:,id]

        # versus eBτ a^2/ħ^2
        if isequal(flag,"versusB")
            # pl1, = plot(abs2.(eBτs[1:length(eBτs)÷2])*1e8,ρ1,"r-^",ms=4)
            pl2, = plot(abs2.(eBτs[1:length(eBτs)÷2])*1e8,ρ2 .- ρ2[end],"-o",ms=4,label="$(νs[iν])")
            xlabel(L"$\left(eBτ\frac{a^2}{ħ^2}\times 10^4\right)^2$")
        elseif isequal(flag,"quadratic")
            # pl1, = plot(abs2.(ωcτ_avg[1:length(eBτs)÷2]),ρ1,"r-^",ms=4)
            pl2, = plot(abs2.(ωcτ_avg[1:length(eBτs)÷2]),(ρ2.- ρ2[end])./(ρ2[end]),"-o",ms=4,label=@sprintf("%.2f",4νs[iν]-4))
            xlabel(L"(ω̄_cτ)^2")
        elseif isequal(flag,"linear")
            # pl1, = plot(abs.(ωcτ_avg[1:length(eBτs)÷2]),ρ1,"r-^",ms=4)
            pl2, = plot(abs.(ωcτ_avg[1:length(eBτs)÷2]),ρ2.- ρ2[end],"-o",ms=4,label="$(νs[iν])")
            xlabel(L"ω̄_cτ")
        end
    end
    
    ylabel(L"$(ρ(H)-ρ(H_0))/ρ(H_0)$")
    # legend([pl1,pl2],[L"$ρ_{1}$",L"$ρ_{2}$"])
    # title("ν=$(ν)")
    legend()
    tight_layout()
    display(fig)
    close(fig)
end

# --------------------------------------------------------------------------------------------------------- #
## plot lattice vector length vs strain angle 
function plot_uc_vector()
    ϕs = collect(0:60)
    a1 = zeros(Float64,length(ϕs))
    a2 = zeros(Float64,length(ϕs))
    a3 = zeros(Float64,length(ϕs))
    area = zeros(Float64,length(ϕs))
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        params = Params(ϵ=0.0068,φ=ϕ*π/180,Da=-4100.0,_hetero=true)
        initParamsWithStrain(params) 
        a1[iϕ] = abs(params.a1)
        a2[iϕ] = abs(params.a2)
        a3[iϕ] = abs(params.a1-params.a2)
        area[iϕ] = imag(params.a1'*params.a2)
    end
    fig = figure(figsize=(2.8,2.8))
    aCC = 0.246
    plot(ϕs,a1*aCC,"r-",ms=4,label=L"L_1")
    plot(ϕs,a2*aCC,"b-",ms=4,label=L"L_2")
    plot(ϕs,a3*aCC,"g-",ms=4,label=L"L_3")
    # plot(ϕs,area*aCC^2,"b-",ms=4)
    print(area*aCC^2)
    # xlim([-2,63])
    # ylim([8,11.5])
    # yticks(8:11)
    # legend()
    xlabel(L"$\varphi\ (^\circ)$")
    # ylabel("Unit Vector Length (nm)")
    ylabel("Moire area")
    tight_layout()
    savefig(joinpath(fpath,"test/moire_area_vs_phi.pdf"),transparent=true)
    display(fig)
    close(fig)
end
# plot_uc_vector()

##
function dirac_splitting()
    ϵDirac = []
    for ϕ in ϕs
        Hk = readdlm(joinpath(fpath,"test/strain_003_phi$(ϕ).txt"))
        δϵk = abs.( Hk[1,:] .- Hk[2,:])
        id = sortperm(δϵk)
        if norm(δϵk[id[1:2]]) < 2
            push!(ϵDirac,abs( (Hk[1,id[1]] + Hk[2,id[1]])/2 - (Hk[1,id[2]] + Hk[2,id[2]])/2 ))
        else
            println("possible accuracy determining Dirac splitting ",ϕ)
            println(norm(δϵk[id[1:2]]))
        end
    end

    fig = figure(figsize=(4,3))
    plot(ϕs,ϵDirac,"b-o",ms=4)
    xlabel(L"$ϕ\ (^\circ)$")
    ylabel("δ (meV)")
    ylim([0,10])
    tight_layout()
    # savefig(joinpath(fpath,"test/dirac_splitting_upper_bound.pdf"),transparent=true)
    display(fig)
    close(fig)
    return nothing

end
# dirac_splitting()

## contour energy maps positive energy 
function plot_contours()
    ϵ2 = reshape(blk.Hk[2,:],Latt.lk,Latt.lk)
    ϵ1 = reshape(blk.Hk[1,:],Latt.lk,Latt.lk)
    kvec = reshape(Latt.kvec,Latt.lk,Latt.lk) ./ abs(8π/sqrt(3)*sin(1.38π/360))
    fig = figure(figsize=(3,2))
    ϵmax = maximum(ϵ2)
    α0, α1,α2 =0.,0.341,0.595
    # pcolor(real(kvec),imag(kvec),ϵ2,cmap="Greens")
    
    ν2 = [4* sum((sign.(α0*ϵmax .-blk.Hk ).+1)./2)/(Latt.lk^2) - 4;
    4* sum((sign.(α1*ϵmax .-blk.Hk ).+1)./2)/(Latt.lk^2) - 4;
    4* sum((sign.(α2*ϵmax .-blk.Hk ).+1)./2)/(Latt.lk^2) - 4; ]
    
    # pl = contour(real(kvec),imag(kvec),ϵ2,levels=[α0*ϵmax,α1*ϵmax, α2*ϵmax])
    # ν2str = Dict(pl.levels[i]=> @sprintf("%1.2f",ν2[i]) for i in eachindex(pl.levels))
    # clabel(pl,pl.levels,fmt=ν2str,inline=true,fontsize=6)
    pl = contour(real(kvec),imag(kvec),ϵ2,levels=[-2,0],colors=["g","r"])
    pl = contour(real(kvec),imag(kvec),ϵ1,levels=[-2,0],colors=["g","r"])
    tight_layout()
    axis("equal")
    title("ϕ=$(ϕ)")
    display(fig)
    
    # savefig("vanhove_phi$(ϕ).pdf",transparent=true)
    close(fig)
end
# plot_contours()

# get van hove points from dos 
function get_dos(ϵ,γ=0.1)
    Evals = range(minimum(ϵ),maximum(ϵ),600)
    dos = zeros(length(Evals))
    for iE in eachindex(Evals)
        dos[iE] =  1/π * sum( γ ./(γ^2 .+ (Evals[iE] .- ϵ).^2 ) ) / length(ϵ)
    end
    return Evals, dos
end

function get_vanhove_from_dos(Hk,Latt)
    ϵ2 = reshape(Hk,:,Latt.lk,Latt.lk)

    Es= [-16.7 21 27.]
    Evals, dos = get_dos(ϵ2,0.4)
    νs = [4* sum((sign.(Es[iE] .-blk.Hk ).+1)./2)/(Latt.lk^2) - 4 for iE in eachindex(Es)]
    fig = figure(figsize=(4,2.6))
    plot(Evals,dos,"k-")
    axvline(Es[1])
    axvline(Es[2])
    axvline(Es[3])
    xlabel("E")
    xlim([-50,50])
    ylabel("DOS")
    yticks([])
    tight_layout()
    # savefig("vanhove_schematic.pdf",transparent=true)
    display(fig)
    close(fig)
    println(Es)
    println(νs)
end
# get_vanhove_from_dos(blk.Hk,Latt)

##
function plot_vanhove(flag="energy")
    df = DataFrame(CSV.File(joinpath(fpath,"test/vanhove_data.csv")))
    fig = figure(figsize=(4,3))
    if isequal(flag,"energy")
        ϵ1 = [df[!,"E1"][1:5];df[!,"E2"][6:21];df[!,"E3"][22:end]]
        ϵ2 = [df[!,"E2"][1:5];df[!,"E1"][6:35];df[!,"E2"][36:end]]
        ϵ3 = [df[!,"E3"][1:21];df[!,"E2"][22:35];df[!,"E1"][36:end]]
        plot(df[!,"Angle"],ϵ1,"r-",ms=4)
        plot(df[!,"Angle"],ϵ2,"b-",ms=4)
        plot(df[!,"Angle"],ϵ3,"g-",ms=4)
        xlabel(L"$ϕ\ (^\circ)$")
        ylabel("E (meV)")
        # axhline(16.5)
        # axvline(60)
        tight_layout()
        savefig("vanhove_energy_vs_phi.pdf",transparent=true)
    elseif isequal(flag,"filling")
        nu1 = [df[!,"nu1"][1:5];df[!,"nu2"][6:21];df[!,"nu3"][22:end]]
        nu2 = [df[!,"nu2"][1:5];df[!,"nu1"][6:35];df[!,"nu2"][36:end]]
        nu3 = [df[!,"nu3"][1:21];df[!,"nu2"][22:35];df[!,"nu1"][36:end]]
        plot(df[!,"Angle"],nu1,"r-",ms=4)
        plot(df[!,"Angle"],nu2,"b-",ms=4)
        plot(df[!,"Angle"],nu3,"g-",ms=4)
        xlabel(L"$ϕ\ (^\circ)$")
        ylabel("ν")
        ylim([0,4])
        # axhline(2.8)
        # axhline(1.35)
        # axvline(60)
        tight_layout()
        # savefig("vanhove_filling_vs_phi.pdf",transparent=true)
    end
    display(fig)
    close(fig)
end

# plot_vanhove("energy")