using Arpack
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
    Hk::Matrix{Float64} # flat band energies nflatxlk^2
    Uk::Matrix{ComplexF64} # gauge-fixed eigenvectors 4lg^2 x nflat lk^2
    Vx::Matrix{Float64} # velocity 
    Vy::Matrix{Float64}
    
    C2T::Matrix{Float64} # Unitary part of the C2T symmetry 
    Ph::Matrix{Float64} # Particle-hole symmetry 

    _σrotation::Bool # flag on whether to add σ matrix rotation 

    Msrk::Array{Float64,2} # nflat x lk^2, Self-rotating magnetization 
    Ωnk::Array{Float64,2} # nflat x lk^2 , Berry curvature, μ part of the chern magnetization
    Mcpartk::Array{Float64,2} # nflat x lk^2, -ϵ_nk part of the chern magnetization

    HBM() = new()
end

@inline function dirac(k::ComplexF64,θ0::Float64) ::Matrix{ComplexF64}
    return  -abs(k)*[0 exp(-1im*(angle(k)-θ0));exp(1im*(angle(k)-θ0)) 0]
end

@inline function diracvalleyKprime(k::ComplexF64,θ0::Float64) ::Matrix{ComplexF64}
    return  -abs(k)*[0 exp(1im*(angle(k)-θ0));exp(-1im*(angle(k)-θ0)) 0]
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
        H[(4idg[ig]-3):(4idg[ig]-2),(4idg[ig]-3):(4idg[ig]-2)] = params.vf*dirac(k1[1]+1im*k1[2],0.0) .- (params.Da * ∇u)*σ0 .+ params.δ*σz
        H[(4idg[ig]-1):(4idg[ig]),(4idg[ig]-1):(4idg[ig])] = params.vf*dirac(k2[1]+1im*k2[2],0.0) .+ (params.Da * ∇u)*σ0 #.+ params.δ*σz
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
    blk.nflat = 30 # take 10 bands above and below
    
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
        idx_mid = (4*blk.lg^2)÷2
        F = eigen(Hermitian(H),(idx_mid-blk.nflat ÷2+1):(idx_mid+blk.nflat÷2))
        vals, vecs = F.values, F.vectors
        # C2T is broken for hBN alignment
        # vecs = vecs + blk.C2T*conj(vecs)
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
    
    ## orbital magnetic moment related
    blk.Msrk = zeros(Float64,blk.nflat,length(Latt.kvec)) # nflat x lk^2, Self-rotating magnetization 
    blk.Ωnk = zeros(Float64,blk.nflat,length(Latt.kvec)) # nflat x lk^2 , Berry curvature, μ part of the chern magnetization
    blk.Mcpartk = zeros(Float64,blk.nflat,length(Latt.kvec)) # nflat x lk^2, -ϵ_nk part of the chern magnetization
    for ik in eachindex(Latt.kvec), iband in 1:blk.nflat 
        _un = view(blk.Uk,:,blk.nflat*(ik-1)+iband)
        _ϵn = blk.Hk[iband,ik]
        for jband in 1:blk.nflat 
            if jband != iband 
                _um = view(blk.Uk,:,blk.nflat*(ik-1)+jband)
                _ϵm = blk.Hk[jband,ik]
                _ker = params.vf^2 * (_un'*σx*_um) * (_um'*σy*_un) / (_ϵn - _ϵm)^2
                blk.Msrk[iband,ik] +=  -imag(_ker * (_ϵn - _ϵm) )
                blk.Ωnk[iband,ik] += -imag(2*_ker)
                blk.Mcpartk[iband,ik] += -imag(-2*_ϵn * _ker)
            end
        end
    end

    return nothing
end

function compute_orbital_magnetization(nband::Int,blk::HBM,area_moire::Float64,area_moire_k::Float64;nμs=200)
    # compute orbtial magnetic moment if chemical potential is inside the nth band 
    ϵn = blk.Hk[nband,:]
    μs = range(0.99*minimum(ϵn),0.99*maximum(ϵn),nμs)
    νs = [sum((sign.(μs[iμ] .- ϵn) .+ 1)./2)/length(ϵn) for iμ in eachindex(μs)]
    msr = zeros(Float64,length(μs))
    mc = zeros(Float64,length(μs))
    mtot = zeros(Float64,length(μs))
    ee = 1.602e-19 * 1e-3 
    hbar = 1.054571817e-34
    me = 9.1093837e-31
    ϵa = hbar^2/(2me*area_moire*(2.46e-10)^2)
    coeff = ee/(ϵa) *area_moire_k /(2π)^2  #*area_moire 
    for iμ in eachindex(μs)
        Θμϵ = (sign.(μs[iμ] .- ϵn) .+ 1)./2
        msr[iμ] = sum(Θμϵ .* view(blk.Msrk,nband,:)) /length(ϵn)
        mc[iμ] = μs[iμ] * sum(Θμϵ .* view(blk.Ωnk,nband,:)) /length(ϵn) + sum(Θμϵ .* view(blk.Mcpartk,nband,:)) /length(ϵn)
        mtot[iμ] = msr[iμ] + mc[iμ]
    end
    return μs,νs,msr*coeff,mc*coeff,mtot*coeff
end