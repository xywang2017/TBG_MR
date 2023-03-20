# plot ρxx, ρyy vs field eBτ
function plot_ρs(transport_coefficients,fname)
    fig = figure(figsize=(5,3))
    pl1, pl2 = 0,0
    νs = [4transport_coefficients[i][1]-4 for i in 1:size(transport_coefficients,1)]
    
    ρxx = [transport_coefficients[i][2][1,1] for i in 1:size(transport_coefficients,1)]
    ρyy = [transport_coefficients[i][2][2,2] for i in 1:size(transport_coefficients,1)]
    
    
    pl1, = plot(νs,ρxx,"-",ms=3)
    pl2, = plot(νs,ρyy,"-",ms=3)
    
    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    for i in 1:length(data)
        axvline(data[i],ls=":",c="gray")
    end
    # 
    ylim([0,40])
    legend([pl1,pl2],[L"ρ_{xx}",L"ρ_{yy}"])
    # legend(loc="upper left")
    xlabel("ν")
    ylabel(L"$ρ/ \left(ρ_Q \frac{Γ}{ϵ_M} \right)$",rotation=90)
    # ylim([0,200])
    # yscale("log")
    # ylim([1.2,8e4])
    tight_layout()
    savefig("unsymmetrized_rho_B0.pdf",transparent=true)
    display(fig)
    close(fig)
    return nothing
end

function plot_ρs_derivative(transport_coefficients,fname)
    fig = figure(figsize=(5,3))
    pl1, pl2 = 0,0
    νs = [4transport_coefficients[i][1]-4 for i in 1:size(transport_coefficients,1)]
    
    ρxx = [transport_coefficients[i][2][1,1] for i in 1:size(transport_coefficients,1)]
    ρyy = [transport_coefficients[i][2][2,2] for i in 1:size(transport_coefficients,1)]
    
    pl1, = plot((νs[1:(end-1)]+νs[2:end])/2,diff(log.(ρxx))./diff(νs),"-",ms=3)
    pl2, = plot((νs[1:(end-1)]+νs[2:end])/2,diff(log.(ρyy))./diff(νs),"-",ms=3)

    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    for i in 1:length(data)
        axvline(data[i],ls=":",c="gray")
    end
    # 
    ylim([-6,6])
    legend([pl1,pl2],[L"ρ_{xx}",L"ρ_{yy}"])
    # legend(loc="upper left")
    xlabel("ν")
    ylabel(L"$d\ln ρ/ d\nu$",rotation=90)
    # ylim([0,200])
    # yscale("log")
    # ylim([1.2,8e4])
    tight_layout()
    savefig("unsymmetrized_rho_B0_derivative.pdf",transparent=true)
    display(fig)
    close(fig)
    return nothing
end

## symmetrization wrt B to -B
function B0symmetrization(transport_coefficients)
    ρ = zeros(2,2,size(transport_coefficients,1))
    νs = zeros(size(transport_coefficients,1))
    for iμ in 1:size(transport_coefficients,1)
        νs[iμ] = transport_coefficients[iμ][1]
        ρ[:,:,iμ] = transport_coefficients[iμ][2] 
    end
    # principle axis
    Uθ = zeros(2,2,size(transport_coefficients,1))
    δθ = zeros(size(transport_coefficients,1))
    for iμ in 1:size(transport_coefficients,1)
        if (abs(ρ[1,2,iμ]-ρ[2,1,iμ])>1e-6)
            println("error with symmetrization ρ[1,2]≂̸ρ[2,1]")
        end
        ρ11,ρ12,ρ22 = ρ[1,1,iμ],ρ[1,2,iμ],ρ[2,2,iμ]
        ω = sqrt((ρ11-ρ22)^2/4+ρ12^2)
        ϵminus = (ρ11+ρ22)/2 - ω
        # δθ[iτ,iμ] = - atan((ρ11-ϵminus)/ρ12)  # (-π/2,π/2)
        δθ[iμ] = acos(-sign(ρ12)*sqrt(0.5*(1-(ρ11-ρ22)/(2ω))))
        Uθ[:,:,iμ] = Float64[cos(δθ[iμ]) -sin(δθ[iμ]);sin(δθ[iμ]) cos(δθ[iμ])]
        ρ[:,:,iμ] = Uθ[:,:,iμ]' * ρ[:,:,iμ] * Uθ[:,:,iμ]
        
        if (abs(ρ[1,2,iμ])+abs(ρ[2,1,iμ])) > 1e-5
            println("error with diagonalization using U")
        end
    end
    return νs,ρ,δθ
end

function plot_symmetrized_transport(νs,ρ,fname)
    fig = figure(figsize=(5,3))
    
    # pl1, = plot(4νs.-4,ρ[1,1,:],"-",ms=3)
    # pl2, = plot(4νs.-4,ρ[2,2,:],"-",ms=3)
    pl2, = plot(4νs.-4,-(1 ./ ρ[2,2,:].- 1 ./ρ[1,1,:]).*ρ[1,1,:],".",ms=3)

    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    for i in 1:length(data)
        axvline(data[i],ls=":",c="gray")
    end

    xlabel("ν")
    # ylabel(L"$ρ/ \left(ρ_Q \frac{Γ}{ϵ_M} \right)$")
    ylabel(L"δσ/σ")
    # legend([pl1,pl2],[L"$ρ_{1}$",L"$ρ_{2}$"],loc="upper left")
    ylim([0,1])
    # yscale("log")
    # xlim([-1,1])
    tight_layout()
    # savefig("symmetrize_rho.pdf",transparent=true)
    display(fig)
    close(fig)
end


##
function plot_principle_strain_axis(νs,δθ,fname)
    fig = figure(figsize=(5,3))
    l0,l1,l2,l3=0,0,0,0

    plot(4νs.-4,δθ*180/π,"-",ms=5)
    
    l1,=plot(4νs.-4,ones(length(νs))*angle(params.a1)*180/π,"--")
    l2,=plot(4νs.-4,ones(length(νs))*angle(params.a2)*180/π,"--")
    l3,=plot(4νs.-4,ones(length(νs))*angle(params.a2-params.a1)*180/π,"--")
    
    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    for i in 1:length(data)
        axvline(data[i],ls=":",c="gray")
    end
    legend([l1,l2,l3],[L"a_1",L"a_2",L"a_3"],loc="upper left")
    xlabel("ν")
    ylabel("δθ")
    yticks([0,60,120],[L"$0$",L"$\frac{π}{3}$",L"$\frac{2π}{3}$"])
    ylim([0,185])
    tight_layout()
    # savefig("principal_axis.pdf",transparent=true)
    display(fig)
    close(fig)
end