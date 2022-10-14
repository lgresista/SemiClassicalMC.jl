struct Configuration
    latticename         :: String                            # e.g. "triangular"
    L                   :: Int64                             # linear extent (in unitcells)
    nSites              :: Int64                             # number of sites
    unitVectors         :: Vector{Vector{Float64}}           # unit vectors of unitcell
    basis               :: Vector{Vector{Float64}}           # Basis sites in unit cell
    project             :: Matrix{Int64}                     # Matrix, so that r(i)-r(j)=getRs(position, basis)[projet[i, j]] (for efficient calculation of correlations)
    positions           :: Vector{Vector{Float64}}           # positions of each site
    siteLabels          :: Vector{Int64}                     # optional site label (for the calculation of sublattice magnetizations)
    interactionSites    :: Vector{Vector{Int64}}             # interactionSites[i] is a vector containing the indices of each site that interacts with site i
    interactionLabels   :: Vector{Vector{Int64}}             # interactionLabels[i] is a vector containing the index of the type of interaction (in interactions)
    interactions        :: Vector{Interaction}               # Stores different types of interactions (nn, nnn, bond depent)
    onsiteInteraction   :: DiagInteraction                   # Onsite interaction (only diagonal to stay Hermitian)
    d                   :: Int64                             # Dimension of local Hilbert space
    generators          :: Vector{Matrix{Complex{Float64}}}  # Basis of generators T of su(4) in matrix form (dxd matrices)
    generatorsSq        :: Vector{Matrix{Complex{Float64}}}  # Square of generators T² (for onsite interaction)
    state               :: Matrix{Complex{Float64}}          # state[:, i] contains local state at site i (complex vector of size d)
    spinExpectation     :: Matrix{Float64}                   # spinExpectation[:, i] constaints <Ψᵢ|T|Ψᵢ> for current state
    spinSqExpectation   :: Matrix{Float64}                   # spinSqExpectation[:, i] constaints <Ψᵢ|T²|Ψᵢ> for current state
end

#Initialization only with necesarry parameters
function Configuration( 
    latticename        :: String,
    L                  :: Int64,
    unitVectors        :: Vector{Vector{Float64}},
    basis              :: Vector{Vector{Float64}},
    positions          :: Vector{Vector{Float64}},
    siteLabels         :: Vector{Int64},
    interactionSites   :: Vector{Vector{Int64}},
    interactionLabels  :: Vector{Vector{Int64}},
    interactions       :: Vector{Interaction},
    onsiteInteraction  :: DiagInteraction,
    generators         :: Vector{Matrix{Complex{Float64}}},
    state              :: Matrix{Complex{Float64}}
)                      :: Configuration

    nSites = length(positions)
    d = size(state, 1)
    project = getProject(positions, L, unitVectors, basis)
    generatorsSq = generators .^2
    spinExpectation = hcat([computeSpinExpectation(state[:, i], generators) for i in 1:nSites]...)
    spinSqExpectation =  hcat([computeSpinExpectation(state[:, i], generatorsSq) for i in 1:nSites]...)

    Configuration(latticename, L, nSites, unitVectors, basis, project, positions, siteLabels, interactionSites,
                  interactionLabels, interactions, onsiteInteraction, d, generators, generatorsSq, 
                  state, spinExpectation, spinSqExpectation)
end

#### Get functions ###

@inline function getLatticename(cfg :: Configuration) :: String
    return cfg.latticename
end

@inline function getL(cfg :: Configuration) :: Int64
    return cfg.l
end

@inline function Base.length(cfg :: Configuration) :: Int64
    return cfg.nSites
end

@inline function getUnitVectors(cfg :: Configuration) :: Vector{Vector{Float64}}
    return cfg.unitVectors
end

@inline function getBasis(cfg :: Configuration) :: Vector{Vector{Float64}}
    return cfg.basis
end

@inline function getProject(cfg :: Configuration) :: Matrix{Int64}
    return cfg.project
end

@inline function getPosition(cfg :: Configuration, i :: Int64) :: Vector{Float64}
    return cfg.positions[i]
end

@inline function getSiteLabel(cfg :: Configuration, i :: Int64) :: Int64
    return cfg.siteLabels[i]
end

@inline function getInteractionSites(cfg :: Configuration, i :: Int64) :: Vector{Int64}
    return cfg.interactionSites[i]
end

@inline function getInteractionLabels(cfg :: Configuration, i :: Int64) :: Vector{Int64}
    return cfg.interactionLabels[i]
end

@inline function getInteraction(cfg :: Configuration, label:: Int64) :: Interaction
    return cfg.interactions[label]
end

@inline function getOnsiteInteraction(cfg :: Configuration) :: DiagInteraction
    return cfg.onsiteInteraction
end

@inline function getD(cfg :: Configuration) :: Int64
    return cfg.d
end

@inline function getGenerators(cfg :: Configuration) :: Vector{Matrix{Complex{Float64}}}
    return cfg.generators
end

@inline function getGeneratorsSq(cfg :: Configuration) :: Vector{Matrix{Complex{Float64}}}
    return cfg.generatorsSq
end

@inline function getState(cfg :: Configuration, i :: Int64) :: Vector{Complex{Float64}}
    return cfg.state[:, i]
end

@inline function getSpinExpectation(cfg :: Configuration) :: Matrix{Float64}
    return cfg.spinExpectation
end

@inline function getSpinExpectation(cfg :: Configuration, i :: Int64) :: Vector{Float64}
    return cfg.spinExpectation[:, i]
end

@inline function getSpinSqExpectation(cfg :: Configuration, i :: Int64) :: Vector{Float64}
    return cfg.spinSqExpectation[:, i]
end

#### Energy and energy difference
function computeSpinExpectation(state :: Vector{Complex{Float64}}, generators :: Vector{Matrix{Complex{Float64}}}) :: Vector{Float64}
    return [real(dot(state, generators[mu], state)) for mu in 1:length(generators)]
end

function getEnergy(cfg :: Configuration) :: Float64
    E = 0.0
    for i in 1:length(cfg)
        Tsq = getSpinSqExpectation(cfg, i)
        E += onsiteEnergy(Tsq, getOnsiteInteraction(cfg))

        interactionSites  = getInteractionSites(cfg, i)
        interactionLabels = getInteractionLabels(cfg, i)
        
        T_i = getSpinExpectation(cfg, i)
        for j in eachindex(interactionSites)
            T_j = getSpinExpectation(cfg, interactionSites[j])
            M = getInteraction(cfg, interactionLabels[j])
            E += exchangeEnergy(T_i, M, T_j) / 2
        end
    end
    return E
end

function getEnergyDifference(cfg :: Configuration, i :: Int64, newT :: Vector{Float64}, newTsq :: Vector{Float64}) :: Float64
    dE = 0.0
    #Onsite Interaction
    Tsq = getSpinSqExpectation(cfg, i)
    dE += onsiteEnergy(newTsq, getOnsiteInteraction(cfg)) - 
          onsiteEnergy(Tsq, getOnsiteInteraction(cfg))

    #Exchange Interaction
    T = getSpinExpectation(cfg, i)
    interactionSites  = getInteractionSites(cfg, i)
    interactionLabels = getInteractionLabels(cfg, i)
    dT = newT .- T

    for j in eachindex(interactionSites)
        T_j = getSpinExpectation(cfg, interactionSites[j])
        M = getInteraction(cfg, interactionLabels[j])
        dE += exchangeEnergy(dT, M, T_j)
    end
    return dE
end

#Different site labels
function getSiteLabels(
    l           :: Lattice{S,B,U},
    latticename :: String
    )           :: Vector{Int64} where {D,N,LB,U,S<:AbstractSite{Int64,D},B<:AbstractBond{LB,N}}
    if latticename == "triangular"
        return getSiteLabelsTriangular(l)
    else
        return ones(Int64, length(l.sites))
    end
end

function getSiteLabelsTriangular(l :: Lattice{S,B,U}
    ) :: Vector{Float64} where {D,N,LB,U,S<:AbstractSite{Int64,D},B<:AbstractBond{LB,N}}
    lattice_vectors = normalize.(l.lattice_vectors)
    
    basis = [lattice_vectors[1] + lattice_vectors[2], 2*lattice_vectors[2] - lattice_vectors[1]]
    siteLabels = zeros(length(l.sites))

    for i in 1:length(l.sites)
        v1 = point(site(l, i))
        v2 = v1 + lattice_vectors[1]

        sol1 = getPrefactors(v1, basis)
        sol2 = getPrefactors(v2, basis)

        if sol1 .% 1 == [0.0, 0.0]
            siteLabels[i] = 1
        elseif sol2 .% 1 == [0.0, 0.0]
            siteLabels[i] = 2
        else
            siteLabels[i] = 3
        end
    end
    return siteLabels
end

# Functions to get project matrix (for efficient all-to-all correlation storage)
function getPrefactors(v :: Vector{<:Number}, 
    basis :: Vector{<:Vector{<:Number}})
    A  = [dot(e1, e2) for e1 in basis, e2 in basis]
    b  = [dot(e, v) for e in basis]
    return round.((A \ b), digits = 2)
end

#Return all possible "distances" between sites in the lattice
function getRs(positions :: Vector{Vector{Float64}}, basis :: Vector{Vector{Float64}}) :: Vector{Vector{Float64}}
    rs = vcat([positions .- Ref(positions[1] .+ b) for b in basis]...)    
end

function getProject(
    positions   :: Vector{Vector{Float64}}, 
    L           :: Int, 
    unitVectors :: Vector{Vector{Float64}}, 
    basis       :: Vector{Vector{Float64}}
    )           :: Matrix{Int64}

    rs = getRs(positions, basis)
    
    project = zeros(Int64, length(positions), length(positions))

    for i in 1:length(positions)
        for j in 1:length(positions)
            r = periodicDistance(positions[i] .- positions[j], unitVectors, basis, L)
            project[i, j] = findfirst(isapprox.(Ref(r), rs; atol = 1e-10))
        end
    end
    return project
end

function periodicDistance(
    r           :: Vector{Float64}, 
    unitVectors :: Vector{Vector{Float64}},
    basis       :: Vector{Vector{Float64}}, 
    L           :: Int
    ) :: Vector{Float64}
    
    bs = [basis[i] - basis[j] for i in eachindex(basis) for j in eachindex(basis) if i!=j]
    bs = [[0.0, 0.0]; bs]
    
    for b in bs
        
        ns = getPrefactors(r .- b, unitVectors)
  
        isint = true
        for n in ns
            if isinteger(n) == false
                isint = false
                break
            end
        end
        if isint == true
            for i in eachindex(ns)
                ns[i] = ((ns[i] %L) + L) %L
            end
            return b .+ sum(ns .* unitVectors)
        end
    end
    @error "r does not seem to lie on lattice"
end