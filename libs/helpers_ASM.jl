## symmetrization wrt B to -B
function Bsymmetrization(transport_coefficients,eBτs,θ0,idx=[1,2])
    ρ = zeros(2,2,1,size(transport_coefficients,1))
    δρ = zeros(2,2,1,size(transport_coefficients,1))
    νs = zeros(size(transport_coefficients,1))
    U = Float64[cos(θ0) -sin(θ0);sin(θ0) cos(θ0)]
    for iμ in 1:size(transport_coefficients,1)
        νs[iμ] = transport_coefficients[iμ][1]
        for iτ in 1:1
            ρ[:,:,iτ,iμ] = U'*(transport_coefficients[iμ][2][:,:,idx[1]] + transport_coefficients[iμ][2][:,:,idx[2]] )/2 * U
            δρ[:,:,iτ,iμ] = U'*(transport_coefficients[iμ][2][:,:,idx[1]] - transport_coefficients[iμ][2][:,:,idx[2]] )/2 * U
        end
    end
    # δθ = zeros(length(eBτs)÷2,size(transport_coefficients,1))
    return νs,ρ,δρ
end

##
function plot_hall(νs,δρs,fname)
    fig = figure(figsize=(5,3))
    pl1, pl2 = 0,0
    area = abs(imag(params.a1'*params.a2)) 
    plot(4νs.-4,(4νs.-4),"--",c="k")
    plot(4νs.-4,4νs,"--",c="k")
    plot(4νs.-4,4νs.-8,"--",c="k")
    idx_select = 1:1
    cnt = 1
    for iτ in idx_select
        # str = @sprintf("%.2f",ω0τs[iτ])
        str = @sprintf("%.2f",0)
        nH = eBτs[iτ] ./ ( δρs[1][2,1,iτ,:] .*( 2π / (params.vf*params.kb) ) ) * area # /(area * 2.46^2*1e-16)
        # nH =  -δρs[1][2,1,iτ,:]
        pl1, = plot(4νs.-4, nH ,"b-",ms=3,label=str,lw=2)
        cnt += 1
    end

    for iτ in idx_select
        # str = @sprintf("%.2f",ω0τs[iτ])
        str = @sprintf("%.3f",0.025)
        nH = eBτs[iτ] ./ ( δρs[2][2,1,iτ,:] .*( 2π / (params.vf*params.kb) ) ) * area # /(area * 2.46^2*1e-16)
        # nH =  -δρs[2][2,1,iτ,:]
        pl1, = plot(4νs.-4, nH ,"r-",ms=3,label=str,lw=2)
        cnt += 1
    end
    
    # axhline(0)
    df = DataFrame( CSV.File(fname) )
    data = sort(df[!,"$(ϕ)_k_nu"])[[2;3;4;9;10;11]]
    # data = sort(df[!,"$(ϕ)_k_nu"])[[2;9]]
    for i in 1:length(data)
        axvline(data[i],ls=":",c="gray")
    end
    axhline(0)
    legend(loc="upper right",title=L"δB/B = ")
    xlabel("ν")
    ylabel(L"\rm n_H\ (moir\'e\ u.c.)")
    # ylabel(L" ρ_{21}")
    ylim([-6,6])
    # xlim([-4,4])
    tight_layout()
    savefig("var0.pdf",transparent=true)
    display(fig)
    close(fig)
end