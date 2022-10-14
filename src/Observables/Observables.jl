abstract type Observables end

#Include different types of observables
include("Observables_generic.jl")
include("Observables_tg_hbn.jl")


#Initialization of observables based on type and configuration
function initializeObservables(obstype :: Type{T}, cfg :: Configuration) :: T where T <: Observables
    obstype(cfg)
end

# Some relevant observables
function getMagnetization(cfg :: Configuration) :: Vector{Float64}
    return dropdims(sum(getSpinExpectation(cfg), dims = 2)/length(cfg), dims = 2)
end

function getSublatticeMagnetization(cfg :: Configuration) :: Matrix{Float64}
    labels = sort(unique(cfg.siteLabels))
    m      = zeros(15, length(labels))

    for (i, label) in enumerate(labels)
        spins = [getSpinExpectation(cfg, j) for j in 1:length(cfg) if getSiteLabel(cfg, j) == label]
        m[:, i] = sum(spins)/length(spins)
    end

    return m
end

function getCorrelations(cfg :: Configuration)
    Ts = getSpinExpectation(cfg)  
    
    
    correlations = zeros(size(Ts, 1), length(getBasis(cfg)) * length(cfg))

    for i in 1:length(cfg)
        for j in 1:length(cfg)
            correlations[:, cfg.project[i,j]] += Ts[:, i] .* Ts[:, j]
        end
    end
    
    return correlations
end

#Scalar chirality in z-direction
function getChirality(cfg :: Configuration, components :: AbstractArray) :: Float64
    κ = 0.0
    for i in 1:length(cfg)
        #Get left and right pointing triangle from site i 
        i2, i3, i4, i5 = getInteractionSites(cfg, i)[[1, 5, 4, 2]]
        
        #Get corresponding spins
        s1 = getSpinExpectation(cfg, i)[components]
        s2 = getSpinExpectation(cfg, i2)[components]
        s3 = getSpinExpectation(cfg, i3)[components]
        s4 = getSpinExpectation(cfg, i4)[components]
        s5 = getSpinExpectation(cfg, i5)[components]

        #compute staggered spin chirality in z direction
        κ += getChirality(s1, s2, s3) - getChirality(s1, s4, s5)
    end
    return κ/length(cfg)/2
end

function getChirality(v1 :: Vector{Float64}, v2 :: Vector{Float64}, v3 :: Vector{Float64}) :: Float64
    return 2/3/sqrt(3) * (v1 × v2 + v2 × v3 + v3 × v1)[3]
end

function getCollinearity(cfg :: Configuration, components :: AbstractArray) :: Float64
    P = 0.0
    N = length(cfg)
    for i in 1:N
        for j in i:N
            P += (dot(getSpinExpectation(cfg, i)[components], getSpinExpectation(cfg, j)[components]))^2
        end
    end
    N_reduced = N/2 * (N+1)
    return 3/2 * (P/N_reduced - 1/3)
end

#Compute structure factor as fouriertransform of sum over components of spin-spin correlations
function computeStructureFactor(correlation, rs, ks; components = eachindex(correlation[:, 1]))
    sf = zeros(length(ks))
    for k in eachindex(ks)
        z = 0.0
        # Compute Fourier transformation at momentum k. 
        for n in 1:length(rs)
            cors = correlation[:, n]
            z += cos(dot(ks[k], rs[n])) * sum(cors[a] for a in components)
        end
        sf[k] = z / length(rs)
    end
    return sf
end

#Compute structure factor for each spin-spin correlation component (T[i]T[i])
function computeStructureFactorVec(correlation, rs, ks)
    sf = zeros(size(correlation, 1), length(ks))
    for k in eachindex(ks)
        z = zeros(size(correlation,1))
        # Compute Fourier transformation at momentum k. 
        for n in 1:length(rs)
            cors = correlation[:, n]
            z .+= cos(dot(ks[k], rs[n])) .* cors
        end
        sf[:, k] .= z ./ length(rs)
    end
    return sf
end

## Momentumspace helper functions
# fold back momentum to first Brillouin zone for given reciprocal lattice vectors
function fold_back!(
    k  :: Vector{Float64},
    k1 :: Vector{Float64},
    k2 :: Vector{Float64}
    )  :: Nothing

    fold    = true
    shifted = deepcopy(k)

    while fold
        fold = false

        for n1 in -1 : 1
            for n2 in -1 : 1
                @. shifted  = k - n1 * k1 - n2 * k2
                abs_shifted = norm(shifted)
                abs_k       = norm(k)

                if abs_shifted < abs_k
                    if abs(abs_shifted - abs_k) > 1e-8
                        k    .= shifted
                        fold  = true 
                    end
                end
            end
        end
    end

    return nothing
end


function getKsInBox(latticename, L; kx_lim = 2π, ky_lim = 2π)
    uc = getUnitcell(Symbol(latticename))
    uc.lattice_vectors .*= L
    ruc = getReciprocalUnitcell(uc)
    rl = getLatticeInBox(ruc, [2 * kx_lim, 2 * ky_lim], zeros(length(ruc.lattice_vectors[1])))
    ks = [site.point for site in rl.sites]
    return ks
end 

function getKsInBZ(lattice, L)
    uc = getUnitcell(:triangular)
    bz = getBrillouinZone(getReciprocalUnitcell(uc))
    corners = bz.corners
    sort!(corners, by = c -> atan(c[1], c[2]))

    lims = maximum(norm.(corners))
    ks = getKsInBox(lattice, L, kx_lim = lims, ky_lim = lims)
    
    nodes = [corners[i][j] for i in 1:length(corners), j in 1:length(corners[1])]
    edges = vcat([[i i+1] for i in 1:length(corners)-1]..., [length(corners) 1])
    verts = [ks[i][j] for i in 1:length(ks), j in 1:length(ks[1])]
    
    inbz = inpoly2(verts, nodes, edges)
    mask = [inbz[i, 1] || inbz[i, 2] for i in 1:length(ks)]
    return ks[mask]
end

