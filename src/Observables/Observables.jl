abstract type Observables end

#Include different types of observables
include("ObservablesGeneric.jl")
include("ObservablesTgHbn.jl")


#Type-stable initialization of observables based on type and configuration
function initializeObservables(obstype :: Type{T}, cfg :: Configuration) :: T where T <: Observables
    obstype(cfg)
end

## Some model-independent observables

#Magnetization over all sites (outputs vector (<T^1>, <T^2>, ...)) for each generator T^i
function getMagnetization(cfg :: Configuration) :: Vector{Float64}
    return sum(getSpinExpectation(cfg))/length(cfg)
end


#Calculate diagonal all-to-all correlations <T^\mu_i><T^\mu_j>, adding all contributions with the same r_i-r_j
@inline function getCorrelations(cfg :: Configuration, project :: Matrix{Int64})

    Ts           = getSpinExpectation(cfg)
    correlations = zeros(length(Ts[1]), length(getBasis(cfg)) * length(cfg))

    for i in 1 : length(cfg)
        Ts_i = Ts[i] 

        for j in 1 : length(cfg)
            Ts_j = Ts[j]
            idx  = project[i, j]
            
            @turbo for k in eachindex(Ts_i)
                correlations[k, idx] += Ts_i[k] * Ts_j[k]
            end
        end
    end
    
    return correlations
end

#Return all possible "r_i - r_j" between sites in the lattice (with periodic boundary conditions)
function getRs(cfg :: Configuration) :: Vector{Vector{Float64}}
    return vcat([cfg.positions .- Ref(cfg.positions[1] .+ b) for b in cfg.basis]...)    
end

#Returns matrix m mapping sites i, j so that r_i - r_j = rs[m[i,j]]
function getProject(cfg :: Configuration) :: Matrix{Int64}

    rs = getRs(cfg)
    
    project = zeros(Int64, length(cfg), length(cfg))

    for i in 1:length(cfg)
        for j in 1:length(cfg)
            r = periodicDistance(i, j, cfg)
            project[i, j] = findfirst(isapprox.(Ref(r), rs; atol = 1e-10))
        end
    end
    return project
end

#Project connecting vector r_i - r_j onto actual lattice sites (beside basis shifts)
function periodicDistance(i :: Int64, j :: Int64, cfg :: Configuration) :: Vector{Float64}
    #Seperate into connecting bravais and basis vector
    db = cfg.basis[cfg.basisLabels[i]] .- cfg.basis[cfg.basisLabels[j]] 
    r = cfg.positions[i] - cfg.positions[j] .- db
    
    #Decompose r = n[1] * unitVectors[1] + n[2] * unitVectors[2]...
    ns = getPrefactors(r, cfg.unitVectors)

    #Fold back onto lattice
    for i in eachindex(ns)
        ns[i] = ((ns[i] % cfg.L) + cfg.L) % cfg.L
    end

    #Add basis difference in the end
    return sum(ns .* cfg.unitVectors) .+ db
end

#Decompose vector v in basis
function getPrefactors(v :: Vector{<:Number}, 
    basis :: Vector{<:Vector{<:Number}})
    A  = [dot(e1, e2) for e1 in basis, e2 in basis]
    b  = [dot(e, v) for e in basis]
    return round.((A \ b), digits = 5)
end

#Get allowed momenta in box around "center" with side lenghts given by dimensions
function getKsInBox(lattice, L, dimensions, center)
    uc = getUnitcell(Symbol(lattice))
    uc.lattice_vectors .*= L
    ruc = getReciprocalUnitcell(uc)
    rl = getLatticeInBox(ruc, dimensions, center)
    ks = [site.point for site in rl.sites]
    return ks
end

#Compute structure factor as fouriertransform of sum over components of spin-spin correlations
function computeStructureFactor(correlation, rs, ks; components = eachindex(correlation[2:end, 1]))
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
