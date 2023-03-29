# using Arpack
using LinearAlgebra
include("Parameters_mod.jl")
include("Lattice_mod.jl")
# --------------------------------------------------------------------------------------------------------------- #
mutable struct HBM
    

    nlocal::Int  # 4 is the local Hilbert space (layerxsublattice)
    lg::Int 
    listG::Matrix{Int}  # indices of the reciprocal lattice vectors g specifying the Hamiltonian 
    gvec::Vector{ComplexF64} # values of gvectors corresponding to the list
    nflat::Int # number of flat bands per valley is 2

    T12::Matrix{ComplexF64} # k-independent part of the Hamiltonian
    Hk::Matrix{Float64} # flat band energies 2xlk^2
    Uk::Matrix{ComplexF64} # gauge-fixed eigenvectors 4lg^2 x nflat lk^2
    Vx::Matrix{Float64} # velocity 
    Vy::Matrix{Float64}
    
    C2T::Matrix{Float64} # Unitary part of the C2T symmetry 
    Ph::Matrix{Float64} # Particle-hole symmetry 

    _σrotation::Bool # flag on whether to add σ matrix rotation 

    HkKbar::Matrix{Float64} # flat band energies in the opposite valley, constructed by symmetry

    HBM() = new()
end

@inline function dirac(k::ComplexF64,θ0::Float64) ::Matrix{ComplexF64}
    return  abs(k)*[0 exp(-1im*(angle(k)-θ0));exp(1im*(angle(k)-θ0)) 0]
end

@inline function diracvalleyKprime(k::ComplexF64,θ0::Float64) ::Matrix{ComplexF64}
    return  -abs(k)*[0 exp(1im*(angle(k)-θ0));exp(-1im*(angle(k)-θ0)) 0]
end

@inline function V(q::ComplexF64) ::Float64
    res = 1e-6
    if abs(q) < res
        Vq = 0
    else
        Vq = 2π/abs(q)
    end
    return Vq
end

function generate_T12(blk::HBM,params::Params)
    # p.b.c. is used 
    idg = reshape(collect(1:blk.lg^2),blk.lg,blk.lg)

    # idg_nn1 = circshift(idg,(-1,0))  # T1 * (|t><b|)
    # idg_nn2 = circshift(idg,(0,-1))  # T2 * (|t><b|)
    # idg_nn12 = circshift(idg,(-1,-1))  # T0 * (|t><b|)

    # per Oskar & Jian choice of g1 and g2
    # idg_nn1 = circshift(idg,(-1,-1))  # T1 * (|t><b|)
    # idg_nn2 = circshift(idg,(0,-1))  # T2 * (|t><b|)
    # idg_nn12 = circshift(idg,(-1,-2))  # T0 * (|t><b|)

    # per Oskar & Jian choice of g1 and g2
    idg_nn1 = circshift(idg,(0,1))  # T1 * (|t><b|)
    idg_nn2 = circshift(idg,(1,1))  # T2 * (|t><b|)
    idg_nn12 = circshift(idg,(0,0))  # T0 * (|t><b|)

    tmp = zeros(ComplexF64,4,blk.lg^2,4,blk.lg^2)

    for ig in eachindex(idg)
        tmp[3:4,idg[ig],1:2,idg_nn1[ig]] = params.T2
        tmp[1:2,idg_nn1[ig],3:4,idg[ig]] = params.T2'

        tmp[3:4,idg[ig],1:2,idg_nn2[ig]] = params.T1
        tmp[1:2,idg_nn2[ig],3:4,idg[ig]] = params.T1'

        tmp[3:4,idg[ig],1:2,idg_nn12[ig]] = params.T0
        tmp[1:2,idg_nn12[ig],3:4,idg[ig]] = params.T0'
    end

    blk.T12 .= reshape(tmp,blk.nlocal*blk.lg^2,blk.nlocal*blk.lg^2)

    return nothing
end


function generate_T12_valleyKprime(blk::HBM,params::Params)
    # p.b.c. is used 
    idg = reshape(collect(1:blk.lg^2),blk.lg,blk.lg)

    # per Oskar & Jian choice of g1 and g2
    idg_nn1 = circshift(idg,(0,-1))  # T1 * (|t><b|)
    idg_nn2 = circshift(idg,(-1,-1))  # T2 * (|t><b|)
    idg_nn12 = circshift(idg,(0,0))  # T0 * (|t><b|)

    tmp = zeros(ComplexF64,4,blk.lg^2,4,blk.lg^2)

    for ig in eachindex(idg)
        tmp[3:4,idg[ig],1:2,idg_nn1[ig]] = params.T2
        tmp[1:2,idg_nn1[ig],3:4,idg[ig]] = params.T2'

        tmp[3:4,idg[ig],1:2,idg_nn2[ig]] = params.T1
        tmp[1:2,idg_nn2[ig],3:4,idg[ig]] = params.T1'

        tmp[3:4,idg[ig],1:2,idg_nn12[ig]] = params.T0
        tmp[1:2,idg_nn12[ig],3:4,idg[ig]] = params.T0'
    end

    blk.T12 .= reshape(tmp,blk.nlocal*blk.lg^2,blk.nlocal*blk.lg^2)

    return nothing
end

function ComputeH(H::Matrix{ComplexF64},blk::HBM,params::Params,k::ComplexF64;ig::Vector{Int}=[0,0])
    """
        Dirac Hamiltonian in the Bloch band basis
    """
    # Note here k only takes values within first mBZ 
    # if k is outside of first BZ, it is labeled by k + ig[1]*blk.g1 + ig[2]*blk.g2
    H .= 0.0 + 0.0im
    σz = ComplexF64[1 0 ; 0 -1]
    σ0 = ComplexF64[1 0; 0 1]
    itr = reshape(collect(1:blk.lg^2),blk.lg,blk.lg)
    idg = view(circshift(itr,(ig[1],ig[2])),:)
    R = params.dθ/2 * Float64[0 -1;1 0]
    ∇u = (params.S[1,1] + params.S[2,2])/2
    # dispersive part
    for ig in 1:blk.lg^2
        qc = blk.gvec[ig]
        kb = k - params.Kb + qc
        kt = k - params.Kt + qc
        if (blk._σrotation==true)
            k1 = (I + R - params.S/2)*[real(kb);imag(kb)]
            k2 = (I - R + params.S/2)*[real(kt);imag(kt)]
            if params._hetero==false
                k1 = (I + R - params.S*params.α)*[real(kb);imag(kb)]
                k2 = (I - R + params.S*(1-params.α))*[real(kt);imag(kt)]
            end
        else
            k1 = [real(kb);imag(kb)]
            k2 = [real(kt);imag(kt)]
        end
        H[(4idg[ig]-3):(4idg[ig]-2),(4idg[ig]-3):(4idg[ig]-2)] = params.vf*dirac(k1[1]+1im*k1[2],0.0) .- (params.Da * ∇u)*σ0
        H[(4idg[ig]-1):(4idg[ig]),(4idg[ig]-1):(4idg[ig])] = params.vf*dirac(k2[1]+1im*k2[2],0.0) .+ (params.Da * ∇u)*σ0
    end
    
    H .= H + blk.T12 - params.μ*I

    return nothing
end

function initHBM(blk::HBM,Latt::Lattice,params::Params;lg::Int=9,_σrotation::Bool=false)
    blk._σrotation = _σrotation
    s0 = Float64[1 0; 0 1]
    s1 = Float64[0 1; 1 0]
    is2 = Float64[0 1; -1 0]

    @assert (lg-1)%2 == 0   # we deal with lg being odd, such that -g and g are related easily
    blk.lg = lg
    blk.listG = zeros(Int,3,lg^2)
    for i2 in 1:lg, i1 in 1:lg
        blk.listG[1,(i2-1)*lg+i1] = i1
        blk.listG[2,(i2-1)*lg+i1] = i2 
        blk.listG[3,(i2-1)*lg+i1] = (i2-1)*lg+i1
    end
    blk.gvec = zeros(ComplexF64,lg^2)
    for ig in 1:lg^2
        blk.gvec[ig] = params.g1 * blk.listG[1,ig] + params.g2 * blk.listG[2,ig]
    end
    G0 = params.g1 * ((lg+1)÷2) + params.g2 * ((lg+1)÷2) # index of 1st Moire BZ is (lg+1)÷2,(lg+1)÷2
    blk.gvec .= blk.gvec .- G0

    # this gives i mu_y I operation in the Bloch basis
    Ig = reverse(Array{Float64}(I,blk.lg^2,blk.lg^2),dims=1)
    blk.Ph = -kron(Ig,kron(is2,s0))

    # this gives C2T eigenstates
    Ig = Array{Float64}(I,blk.lg^2,blk.lg^2)
    blk.C2T = kron(Ig,kron(s0,s1)) # × conj(...)

    blk.nlocal = 4
    blk.nflat = 2  ### !!
    
    blk.Uk = zeros(ComplexF64,blk.nlocal*blk.lg^2,blk.nflat*length(Latt.kvec))
    blk.Hk =zeros(Float64,blk.nflat,length(Latt.kvec))

    blk.T12 = zeros(ComplexF64,blk.nlocal*blk.lg^2,blk.nlocal*blk.lg^2)
    generate_T12(blk,params)

    # temporary container H for each k
    H = zeros(ComplexF64,blk.nlocal*blk.lg^2,blk.nlocal*blk.lg^2)

    for ik in eachindex(Latt.kvec)
        kval = Latt.kvec[ik] #+ 0.5*(1+1im)/Latt.lk
        ComputeH(H,blk,params,kval)
        # Find the smallest eigenvalue and eigenvectors close to zero
        # vals, vecs = eigs(Hermitian(H),nev=blk.nflat,which=:SM)
        flat_band_index = (size(H,2)÷2):(size(H,2)÷2+1)
        vals, vecs = eigen(Hermitian(H),flat_band_index)
        # C2T is broken for hBN alignment
        vecs = vecs + blk.C2T*conj(vecs)
        for i in 1:blk.nflat
            tmp = view(vecs,:,i)
            normalize!(tmp)
        end

        if (norm(imag(vals))<1e-6)
            perm = sortperm(real(vals[:]))
            blk.Uk[:,(blk.nflat*(ik-1)+1):(blk.nflat*ik)] = view(vecs,:,perm)
            blk.Hk[:,ik] = real(vals[perm])
        else
            print("Error with Hermiticity of Hamiltonian!\n")
        end

    end

    # compute velocities 
    Ig = Array{ComplexF64}(I,blk.lg^2,blk.lg^2)
    σx = zeros(ComplexF64,4blk.lg^2,4blk.lg^2)
    σy = zeros(ComplexF64,4blk.lg^2,4blk.lg^2)
    if (blk._σrotation==true)
        R = ComplexF64[0 -params.dθ/2;params.dθ/2 0]
        M1 = I - R - params.S/2 
        M2 = I + R + params.S/2 
        if params._hetero==false
            M1 = I - R - params.S *params.α
            M2 = I + R + params.S*(1-params.α)
        end
        σx1 = M1[1,1] * ComplexF64[0 1;1 0] + M1[1,2] * ComplexF64[0 -1im; 1im 0]
        σy1 = M1[2,1] * ComplexF64[0 1;1 0] + M1[2,2] * ComplexF64[0 -1im; 1im 0]
        σx2 = M2[1,1] * ComplexF64[0 1;1 0] + M2[1,2] * ComplexF64[0 -1im; 1im 0]
        σy2 = M2[2,1] * ComplexF64[0 1;1 0] + M2[2,2] * ComplexF64[0 -1im; 1im 0]
        σx .= kron(Ig,[σx1 0*σx1;0*σx2 σx2])
        σy .= kron(Ig,[σy1 0*σy1;0*σy2 σy2])
    else
        σx .= kron( Ig, kron(ComplexF64[1 0;0 1],ComplexF64[0 1;1 0]) )
        σy .= kron( Ig, kron(ComplexF64[1 0;0 1],ComplexF64[0 -1im;1im 0]) )
    end

    blk.Vx =zeros(Float64,blk.nflat,length(Latt.kvec))
    blk.Vy =zeros(Float64,blk.nflat,length(Latt.kvec))
    for ik in eachindex(Latt.kvec), iband in 1:blk.nflat
        uvec = view(blk.Uk,:,blk.nflat*(ik-1)+iband)
        blk.Vx[iband,ik] = params.vf * real( uvec' * σx * uvec )
        blk.Vy[iband,ik] = params.vf * real( uvec' * σy * uvec )
    end
    
    # blk.HkKbar = zeros(Float64,blk.nflat,length(Latt.kvec))
    # blk.HkKbar .= - circshift(blk.Hk,(1,0))
    return nothing
end

function find_special_points_BM(blk::HBM,Latt::Lattice,params::Params)
    # for each band, there are at a minimum 3 van Hove points, 2 Dirac points, and a band minima/maxima 
    special_points = []
    kvec = reshape(Latt.kvec,Latt.lk,Latt.lk)
    k12grid = reshape(Latt.k1,:,1) .+ 1im*reshape(Latt.k2,1,:)
    for iband in 1:2
        special_points_n = ComplexF64[]
        # fig = figure(figsize=(3,3))
        pl = contour(real(k12grid),imag(k12grid),reshape(blk.Vx[iband,:],Latt.lk,Latt.lk),levels=[0],colors="r")
        # pl0 = contour(real(k12grid),imag(k12grid),reshape(blk.Vy[iband,:],Latt.lk,Latt.lk),levels=[0],colors="b")
        # tight_layout()
        # display(fig)
        # close(fig)

        segs = [pl.allsegs[i] for i in 1:size(pl.allsegs,2)] 
        # matching boundary points
        paths = zeros(ComplexF64,2*size(segs,1)) 
        for i in 1:size(segs,1)
            paths[2i-1] = 2π*[1; -1im]'*segs[i][1,:] # begin k point, in exp(i2π k1 + 2πk2) format
            paths[2i] = 2π* [1; -1im]'*segs[i][end,:] # end k point
        end

        ## find matching begin/end points
        pairs = []
        for i in 1:length(paths)
            dist_x = abs.(exp(1im * real(paths[i])).- exp.(1im *real(paths)) )
            dist_y = abs.(exp(1im * imag(paths[i])).- exp.(1im *imag(paths)) )
            dist = sqrt.(dist_x.^2 .+ dist_y.^2)
            dist[i] = 1e10
            k = argmin(dist)
            if k <= i
                # always try push end of a chain; pair index stored in [seg_no;end; seg_no; begin/end]
                # if [seg1,1,seg2,1] occurs, then push 1, and later pay attention to store the reverse seg 
                if (i-1)%2+1 == 2  
                    push!(pairs,[(i-1)÷2+1;(i-1)%2+1;(k-1)÷2+1;(k-1)%2+1])
                else 
                    push!(pairs,[(k-1)÷2+1;(k-1)%2+1;(i-1)÷2+1;(i-1)%2+1])
                end
            end
        end
        
        # if iband == 2
        #     println(pairs)
        #     fig = figure(figsize=(3,3))
        #     for i in 1:size(segs,1)
        #         plot(segs[i][:,1],segs[i][:,2],"k--",label="$(i)")
        #     end
        #     colors = ["r","b","g","m","y"]
        #     for i in 1:size(pairs,1)
        #         seg1 = pairs[i][1]
        #         seg1_point = pairs[i][2]
        #         seg2 = pairs[i][3]
        #         seg2_point = pairs[i][4]
        #         if seg1_point == 2 && seg2_point == 1
        #             plot([segs[seg1][end,1],segs[seg2][1,1]],[segs[seg1][end,2],segs[seg2][1,2]],"o",c=colors[i])
        #         end
        #     end
        #     legend()
        #     axis("equal")
        #     tight_layout()
        #     display(fig)
        #     close(fig)
        # end

        # println(pairs)
        ## combine paths by gluing matching pairs of end points together 
        path = [] 
        while (! isempty(pairs) )
            chain = []
            i_pair = 1
            ii, id_ii, oo, id_oo = pairs[i_pair]
            if oo == ii # self-closed loop 
                push!(chain,pairs[i_pair])
                deleteat!(pairs,i_pair)
            else # open ended loop needs to be glued via zone boundaries 
                push!(chain,pairs[i_pair])
                deleteat!(pairs,i_pair)
                j = 1
                while (! isempty(pairs) )
                    if pairs[j][1] == oo 
                        push!(chain,pairs[j])
                        _, _, oo, id_oo = pairs[j]
                        deleteat!(pairs,j)
                        j = 1
                        if (oo == ii) # close loop found 
                            break 
                        end
                    else 
                        j = j+1
                    end
                end
            end
            push!(path,chain)
        end
        # println("Found $(size(path,1)) FS segments stored in path")

        # combine glued data into distinct FSs 
        FS = []
        for iseg in 1:size(path,1)
            chain = path[iseg]
            data = Float64[0.0 0.0]
            for j in 1:size(chain,1)
                if chain[j][2] == 2
                    data = [data; segs[chain[j][1]]]
                else
                    data = [data; reverse(segs[chain[j][1]],dims=1)]
                end
            end
            push!(FS,data[2:end,:])
        end
        println("number of disconnected segments found is: ", size(FS,1))

        for iseg in 1:size(FS,1)
            kseg = FS[iseg][:,1]*params.g1 + FS[iseg][:,2]*params.g2
            vx, vy = vk(kseg,params,iband)
            println(norm(vx)/length(vx))  # check how far away vx is from 0 
            sign_vy = sign.(vy)
            for j in 1:length(sign_vy)
                if sign_vy[mod(j,length(sign_vy))+1] == - sign_vy[j]
                    push!(special_points_n,(kseg[j]+kseg[mod(j,length(sign_vy))+1])/2)
                end
            end
        end
        push!(special_points,special_points_n)
    end

    special_point_energies = [] 
    for iband in 1:2 
        special_point_energies_n = Float64[]
        for j in eachindex(special_points[iband])
            kvals = [special_points[iband][j]]
            latt = initLatticeWithKvec(kvals)
            blk = HBM()
            initHBM(blk,latt,params;lg=9)   
            push!(special_point_energies_n,blk.Hk[iband,1])
        end
        push!(special_point_energies,special_point_energies_n)
    end
    return special_points, special_point_energies
end

function vk(kvec::Vector{ComplexF64},params::Params,iband::Int)
    ## init lattice 
    latt = initLatticeWithKvec(kvec)
    blk = HBM()
    initHBM(blk,latt,params;lg=9)
    # compute velocities 
    Ig = Array{ComplexF64}(I,blk.lg^2,blk.lg^2)
    σx = zeros(ComplexF64,4blk.lg^2,4blk.lg^2)
    σy = zeros(ComplexF64,4blk.lg^2,4blk.lg^2)
    if (blk._σrotation==true)
        R = ComplexF64[0 -params.dθ/2;params.dθ/2 0]
        M1 = I - R - params.S/2 
        M2 = I + R + params.S/2 
        if params._hetero==false
            M1 = I - R - params.S*params.α
            M2 = I + R + params.S*(1-params.α)
        end
        σx1 = M1[1,1] * ComplexF64[0 1;1 0] + M1[1,2] * ComplexF64[0 -1im; 1im 0]
        σy1 = M1[2,1] * ComplexF64[0 1;1 0] + M1[2,2] * ComplexF64[0 -1im; 1im 0]
        σx2 = M2[1,1] * ComplexF64[0 1;1 0] + M2[1,2] * ComplexF64[0 -1im; 1im 0]
        σy2 = M2[2,1] * ComplexF64[0 1;1 0] + M2[2,2] * ComplexF64[0 -1im; 1im 0]
        σx .= kron(Ig,[σx1 0*σx1;0*σx2 σx2])
        σy .= kron(Ig,[σy1 0*σy1;0*σy2 σy2])
    else
        σx .= kron( Ig, kron(ComplexF64[1 0;0 1],ComplexF64[0 1;1 0]) )
        σy .= kron( Ig, kron(ComplexF64[1 0;0 1],ComplexF64[0 -1im;1im 0]) )
    end

    vx = zeros(Float64,length(kvec))
    vy = zeros(Float64,length(kvec))
    for ik in eachindex(kvec)
        uvec = view(blk.Uk,:,blk.nflat*(ik-1)+iband)
        vx[ik] = params.vf * real( uvec' * σx * uvec )
        vy[ik] = params.vf * real( uvec' * σy * uvec )
    end
    return vx,vy
end


function generate_contours_torus(μ::Float64,ϵ::Matrix{Float64},kvec::Matrix{ComplexF64})
    # ϵ[kx,ky], kvec[1,2] = k1[i] + 1im k2[j]; periodicity under k1/k2 -> k1/k2 + 1
    # FS is a vector containing all segments of the Fermi surface 
    
    # all segments of contours generated from matplotlib.pyplot.contour command
    
    pl = contour(real(kvec),imag(kvec),ϵ,levels=[μ])
    segs = [pl.allsegs[i] for i in 1:size(pl.allsegs,2)] 

    # colors = ["r","b","g","m"]
    # fig = figure(figsize=(3,3))
    # for i in 1:size(segs,1)
    #     plot(segs[i][:,1],segs[i][:,2],"-",c=colors[i],label="$(i)")
    #     plot(segs[i][1,1],segs[i][1,2],"x",c="k")
    #     plot(segs[i][end,1],segs[i][end,2],"o",c="k")
    # end
    # legend()
    # axis("equal")
    # tight_layout()
    # display(fig)
    # close(fig)

    # matching boundary points
    paths = zeros(ComplexF64,2*size(segs,1)) 
    for i in 1:size(segs,1)
        paths[2i-1] = 2π*[1; -1im]'*segs[i][1,:] # begin k point, in exp(i2π k1 + 2πk2) format
        paths[2i] = 2π* [1; -1im]'*segs[i][end,:] # end k point
    end

    ## find matching begin/end points
    pairs = []
    for i in 1:length(paths)
        dist_x = abs.(exp(1im * real(paths[i])).- exp.(1im *real(paths)) )
        dist_y = abs.(exp(1im * imag(paths[i])).- exp.(1im *imag(paths)) )
        dist = sqrt.(dist_x.^2 .+ dist_y.^2)
        dist[i] = 1e10
        k = argmin(dist)
        if k <= i
            # always try push end of a chain; pair index stored in [seg_no;end; seg_no; begin/end]
            # if [seg1,1,seg2,1] occurs, then push 1, and later pay attention to store the reverse seg 
            if (i-1)%2+1 == 2  
                push!(pairs,[(i-1)÷2+1;(i-1)%2+1;(k-1)÷2+1;(k-1)%2+1])
            else 
                push!(pairs,[(k-1)÷2+1;(k-1)%2+1;(i-1)÷2+1;(i-1)%2+1])
            end
        end
    end
    # println(pairs)
    
    # fig = figure(figsize=(3,3))
    # for i in 1:size(segs,1)
    #     plot(segs[i][:,1],segs[i][:,2],"k--",label="$(i)")
    # end
    # colors = ["r","b","g","m"]
    # for i in 1:size(pairs,1)
    #     seg1 = pairs[i][1]
    #     seg1_point = pairs[i][2]
    #     seg2 = pairs[i][3]
    #     seg2_point = pairs[i][4]
    #     if seg1_point == 2 && seg2_point == 1
    #         plot([segs[seg1][end,1],segs[seg2][1,1]],[segs[seg1][end,2],segs[seg2][1,2]],"o",c=colors[i])
    #     end
    # end
    # legend()
    # axis("equal")
    # tight_layout()
    # display(fig)
    # close(fig)

    # println(pairs)
    ## combine paths by gluing matching pairs of end points together 
    path = [] 
    while (! isempty(pairs) )
        chain = []
        i_pair = 1
        ii, id_ii, oo, id_oo = pairs[i_pair]
        if oo == ii # self-closed loop 
            push!(chain,pairs[i_pair])
            deleteat!(pairs,i_pair)
        else # open ended loop needs to be glued via zone boundaries 
            push!(chain,pairs[i_pair])
            deleteat!(pairs,i_pair)
            j = 1
            while (! isempty(pairs) )
                if pairs[j][1] == oo 
                    push!(chain,pairs[j])
                    _, _, oo, id_oo = pairs[j]
                    deleteat!(pairs,j)
                    j = 1
                    if (oo == ii) # close loop found 
                        break 
                    end
                else 
                    j = j+1
                end
            end
        end
        push!(path,chain)
    end
    # println("Found $(size(path,1)) FS segments stored in path")

    # combine glued data into distinct FSs 
    FS = []
    for iseg in 1:size(path,1)
        chain = path[iseg]
        data = Float64[0.0 0.0]
        for j in 1:size(chain,1)
            if chain[j][2] == 2
                data = [data; segs[chain[j][1]]]
            else
                data = [data; reverse(segs[chain[j][1]],dims=1)]
            end
        end
        push!(FS,data[2:end,:])
    end

    # fine tune FS to remove small differences in k points 
    for iFS in 1:size(FS,1)
        dk = sqrt.(diff(FS[iFS][:,1]).^2 + diff(FS[iFS][:,2]).^2) 
        dkavg = sum(dk)/length(dk)
        points_to_keep = [1; (2:length(dk)+1)[dk .> 0.2dkavg] ]
        FS[iFS] = FS[iFS][points_to_keep,:]
    end
    return FS
end

function boltzmann_characteristics(FS::Vector{Any},params::Params,iband::Int)
    # solve equation for s based on dkx/ds = qvyB and dky/ds = - qvxB
    ## Fermi velocities 
    vx, vy, s = [], [], []
    FSvec = []
    for i in 1:size(FS,1)
        FS_kchain = FS[i][:,1]*params.g1 + FS[i][:,2]*params.g2
        v1, v2 =  vk(FS_kchain,params,iband)
        push!(vx,v1)
        push!(vy,v2)
        push!(FSvec,FS_kchain)
        tmps = zeros(Float64,size(FS[i],1))
        for j in 2:length(FS_kchain)
            δk1_list = (FS[i][j,1]-FS[i][j-1,1]).+[0;-1;1]
            δk1 = δk1_list[argmin(abs.(δk1_list))] 
            δk2_list = (FS[i][j,2]-FS[i][j-1,2]).+[0;-1;1]
            δk2 = δk2_list[argmin(abs.(δk2_list))] 
            δk = δk1 * params.g1 + δk2 * params.g2
            Δsx = real(δk) / (0.5*vy[i][j] + 0.5*vy[i][j-1])
            Δsy = - imag(δk) / (0.5*vx[i][j] + 0.5*vx[i][j-1])
            # if (sign(Δsx)!=sign(Δsy))
            #     println("error with sign of Δs")
            # end
            
            # avoid dividing over extremely small Fermi velocities 
            # if (abs(Δsx)<abs(Δsy))
            #     tmps[j] = tmps[j-1] + Δsx 
            # else
            #     tmps[j] = tmps[j-1] + Δsy
            # end
            signΔs = abs(Δsx) < abs(Δsy) ? sign(Δsx) : sign(Δsy)
            # Δs = signΔs * abs(δk) / sqrt((vx[i][j] + vx[i][j-1])^2/4 + (vy[i][j] + vy[i][j-1])^2/4)
            Δs = signΔs * abs(δk) / (0.5*sqrt(vx[i][j]^2 + vy[i][j]^2) +0.5*sqrt(vx[i][j-1]^2 + vy[i][j-1]^2))
            tmps[j] = tmps[j-1] + Δs
        end
        push!(s,tmps)
    end
    # save("closedFS.jld","vx",vx,"vy",vy,"FS",FSvec,"s",s)
    # fourier transform 
    vxn,vyn = [], []
    ns = collect(-20:20) 
    for i in 1:size(s,1)
        svec = (s[i][2:end] .+ s[i][1:(end-1)] ) ./2
        vxvec = 0.5*(vx[i][2:end].+vx[i][1:(end-1)])
        vyvec = 0.5*(vy[i][2:end].+vy[i][1:(end-1)])
        dsvec = s[i][2:end] .- s[i][1:(end-1)] 
        push!(vxn, [ sum( dsvec .* vxvec .* exp.(-1im * 2π * n * svec/s[i][end]) )/s[i][end] for n in ns])
        push!(vyn, [ sum( dsvec .* vyvec .* exp.(-1im * 2π * n * svec/s[i][end]) )/s[i][end] for n in ns])
    end
    return vxn,vyn,s  # (in units of 1/qB)
end

function computeTransport(vxn::Vector{Any},vyn::Vector{Any},s::Vector{Any},eBτs::Vector{Float64},params::Params)
    # this is transport by summing up contours of a given energy 
    ns = collect(-20:20)
    ωcτs = Float64[]
    # println(2π / s[1][end] * eBτs[1])
    σ_band = zeros(Float64,2,2,length(eBτs))
    σ = Float64[0 0;0 0] 
    for II in eachindex(eBτs)
        eBτ = eBτs[II]
        σ .= 0.0
        tmp = []
        for i in 1:size(s,1) # σ in units of e^2/ħ * τ/ħ
            ωcτ = 2π / s[i][end] * eBτ
            push!(ωcτs,ωcτ)
            σ[1,1] += real( sum(reverse(vxn[i]).*vxn[i] ./ (1 .+ 1im * ns * ωcτ)) ) * abs(s[i][end])/(2π)^2 #* sign(eBτ)
            σ[1,2] += real( sum(reverse(vxn[i]).*vyn[i] ./ (1 .+ 1im * ns * ωcτ)) ) * abs(s[i][end])/(2π)^2 #* sign(eBτ)
            σ[2,1] += real( sum(reverse(vyn[i]).*vxn[i] ./ (1 .+ 1im * ns * ωcτ)) ) * abs(s[i][end])/(2π)^2 #* sign(eBτ)
            σ[2,2] += real( sum(reverse(vyn[i]).*vyn[i] ./ (1 .+ 1im * ns * ωcτ)) ) * abs(s[i][end])/(2π)^2 #* sign(eBτ)
        end
        # in units of ρ_Q * γ/ϵ_M, ϵ_M = v_F k_θ
        # 4 comes from valley and spin if C2 is preserved
        σ_band[:,:,II] = (2π) /( params.vf * params.kb ) * (4*σ)
        # ρ[:,:,II] = inv(4*σ) / (2π) * params.vf * params.kb 
        # nH[II] = eBτ / real(ρ[2,1]) # only works if principle axis
    end
    return σ_band,ωcτs
end

function mainTransport(blk::HBM,Latt::Lattice,params::Params;eBτs::Vector{Float64}=collect(-3:3),nμs::Int=10)
    # k1 + 1im k2
    k12grid = reshape(Latt.k1,:,1) .+ 1im * reshape(Latt.k2,1,:)
    # then define μs and νs
    # W = (maximum(blk.Hk) - minimum(blk.Hk))
    W_min = minimum(blk.Hk)
    W_max = maximum(blk.Hk)
    μstmp = range(0.96*W_min,0.96*W_max,length=10nμs)
    # μs = ( minimum(blk.Hk) .+ W/ nμs * collect(0.5:(nμs-0.5)))
    νstmp = [sum( (sign.(μstmp[i] .- reshape(blk.Hk,2,Latt.lk,Latt.lk)[:,1:(end-1),1:(end-1)]) .+1) ./2)/(Latt.lk-1)^2 for i in eachindex(μstmp)]
    νs = Float64[νstmp[1]]
    μs = Float64[μstmp[1]]
    for i in 2:length(νstmp)
        if (νstmp[i]-νs[end]) >0.002
            push!(νs,νstmp[i])
            push!(μs,μstmp[i])
        end
    end
    # for each μ, and band index, call FS = generate_contours_torus()
    transport_coefficients = []
    σ0 = zeros(Float64,2,2,length(eBτs))
    FS_collect = []
    for iμ in eachindex(μs)
        ωcτs = []
        if abs(4νs[iμ]-4) < 1e-5 || abs(4νs[iμ]-4) > 0.8
            continue 
        end
        println(4νs[iμ]-4)
        # if (4νs[iμ]-4) < 2.3
        #     continue
        # end
        ρ0 = zeros(Float64,2,2,length(eBτs))
        σ0 .= 0.0
        tmpFS = []
        for iband in 1:2 # need to add valley later
            if μs[iμ] < maximum(blk.Hk[iband,:]) && μs[iμ] > minimum(blk.Hk[iband,:])
                FS = generate_contours_torus(μs[iμ],reshape(blk.Hk[iband,:],Latt.lk,Latt.lk), k12grid)
                vxn, vyn, s = boltzmann_characteristics(FS,params,iband)
                σσ,ωcτ = computeTransport(vxn,vyn,s,eBτs,params)
                σ0 .= σ0 + σσ
                push!(ωcτs,ωcτ)
                push!(tmpFS,copy(FS))
            end
        end
        push!(FS_collect,[4νs[iμ]-4,μs[iμ],tmpFS])
        for iτ in eachindex(eBτs)
            ρ0[:,:,iτ] = inv(σ0[:,:,iτ])
        end
        push!(transport_coefficients,[νs[iμ],ρ0,ωcτs])
    end
    return transport_coefficients , FS_collect
end

# ------------------------------- B=0 transport properties ----------------------------- #
function boltzmann_characteristicsB0(FS::Vector{Any},params::Params,iband::Int)
    # solve equation for s based on dkx/ds = qvyB and dky/ds = - qvxB
    ## Fermi velocities 
    vx, vy, s = [], [], []
    for i in 1:size(FS,1)
        FS_kchain = FS[i][:,1]*params.g1 + FS[i][:,2]*params.g2
        v1, v2 =  vk(FS_kchain,params,iband)
        push!(vx,v1)
        push!(vy,v2)
        Tmps = zeros(Float64,length(FS_kchain)-1)
        for j in 2:length(FS_kchain)
            δk1_list = (FS[i][j,1]-FS[i][j-1,1]).+[0;-1;1]
            δk1 = δk1_list[argmin(abs.(δk1_list))] 
            δk2_list = (FS[i][j,2]-FS[i][j-1,2]).+[0;-1;1]
            δk2 = δk2_list[argmin(abs.(δk2_list))] 
            δk = δk1 * params.g1 + δk2 * params.g2
            Tmps[j-1] = abs(δk)
        end
        push!(s,Tmps)
    end
    return vx,vy,s  # vx, vy, and dk parallel 
end

function computeTransportB0(vx::Vector{Any},vy::Vector{Any},s::Vector{Any},params::Params)
    # this is transport by summing up contours of a given energy 
    # println(2π / s[1][end] * eBτs[1])
    σ_band = zeros(Float64,2,2)
    σ = Float64[0 0;0 0] 
    for i in 1:size(s,1)
        vxavg = (vx[i][2:end].+vx[i][1:(end-1)])./2
        vyavg = (vy[i][2:end].+vy[i][1:(end-1)])./2
        vavg = sqrt.(vxavg.^2 .+ vyavg.^2)
        σ[1,1] += sum(s[i] .* vxavg.^2 ./ vavg)  /(2π)^2
        σ[2,2] += sum(s[i] .* vyavg.^2 ./ vavg) /(2π)^2
        σ[1,2] += sum(s[i] .* vxavg.*vyavg./vavg) /(2π)^2
    end
    σ[2,1] = σ[1,2]
    # in units of ρ_Q * γ/ϵ_M, ϵ_M = v_F k_θ
    # 4 comes from valley and spin if C2 is preserved
    σ_band = (2π) /( params.vf * params.kb ) * (4*σ)

    return σ_band
end

function mainTransportB0(blk::HBM,Latt::Lattice,params::Params;nμs::Int=10)
    # k1 + 1im k2
    k12grid = reshape(Latt.k1,:,1) .+ 1im * reshape(Latt.k2,1,:)
    # then define μs and νs
    # W = (maximum(blk.Hk) - minimum(blk.Hk))
    W_min = minimum(blk.Hk)
    W_max = maximum(blk.Hk)
    μstmp = range(0.96*W_min,0.96*W_max,length=50nμs)
    # μs = ( minimum(blk.Hk) .+ W/ nμs * collect(0.5:(nμs-0.5)))
    νstmp = [sum( (sign.(μstmp[i] .- reshape(blk.Hk,2,Latt.lk,Latt.lk)[:,1:(end-1),1:(end-1)]) .+1) ./2)/(Latt.lk-1)^2 for i in eachindex(μstmp)]
    νs = Float64[νstmp[1]]
    μs = Float64[μstmp[1]]
    for i in 2:length(νstmp)
        if abs(νs[end])>0.4 
            if νstmp[i]-νs[end] >0.01
                push!(νs,νstmp[i])
                push!(μs,μstmp[i])
            end
        else  
            if νstmp[i]-νs[end] >0.004
                push!(νs,νstmp[i])
                push!(μs,μstmp[i])
            end
        end
    end
    # for each μ, and band index, call FS = generate_contours_torus()
    transport_coefficients = []
    σ0 = zeros(Float64,2,2)
    for iμ in eachindex(μs)
        println(4νs[iμ]-4)
        if abs(4νs[iμ]-4) < 1e-5 
            continue 
        end
        ρ0 = zeros(Float64,2,2)
        σ0 .= 0.0
        for iband in 1:2 # need to add valley later
            if μs[iμ] < maximum(blk.Hk[iband,:]) && μs[iμ] > minimum(blk.Hk[iband,:])
                FS = generate_contours_torus(μs[iμ],reshape(blk.Hk[iband,:],Latt.lk,Latt.lk), k12grid)
                vx, vy, s = boltzmann_characteristicsB0(FS,params,iband)
                σσ = computeTransportB0(vx,vy,s,params)
                σ0 .= σ0 + σσ
            end
        end
        ρ0[:,:] = inv(σ0)
        push!(transport_coefficients,[νs[iμ],ρ0])
    end
    return transport_coefficients
end

