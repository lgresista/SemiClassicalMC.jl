# Configuration type
struct Configuration
    latticename         :: String                           # e.g. "triangular"
    L                   :: Int64                            # linear extent (in unitcells)
    nSites              :: Int64                            # number of sites
    unitVectors         :: Vector{Vector{Float64}}          # bravais Vectors
    basis               :: Vector{Vector{Float64}}          # basis vectors
    positions           :: Vector{Vector{Float64}}          # positions of each site
    basisLabels         :: Vector{Int64}                    # basis sublattice label of each site
    interactionSites    :: Vector{Vector{Int64}}            # interactionSites[i] is a vector containing the indices of each site that interacts with site i
    interactionLabels   :: Vector{Vector{Int64}}            # interactionLabels[i] is a vector containing the index of the type of interaction (in interactions)
    interactions        :: Vector{SMatrix{16, 16, Float64, 256}} # Stores different types of interactions (nn, nnn, bond depent)
    onsiteInteraction   :: SVector{16, Float64}             # Onsite interaction (only diagonal to stay Hermitian)
    d                   :: Int64                            # Dimension of local Hilbert space
    generators          :: Vector{Matrix{Complex{Float64}}} # Basis of generators T of su(4) in matrix form (dxd matrices)
    generatorsSq        :: Vector{Matrix{Complex{Float64}}} # Square of generators T² (for onsite interaction)
    state               :: Vector{Vector{Complex{Float64}}} # state[:, i] contains local state at site i (complex vector of size d)
    spinExpectation     :: Vector{MVector{16, Float64}}     # spinExpectation[i] constaints <Ψᵢ|T|Ψᵢ> for current state
    spinSqExpectation   :: Vector{MVector{16, Float64}}     # spinSqExpectation[i] constaints <Ψᵢ|T²|Ψᵢ> for current state
    dT                  :: MVector{16, Float64}             # Buffer to hold dT = <T_i_new> - <T_i_old> in local update
    newT                :: MVector{16, Float64}             # Buffer to hold <T_i_new>
    newTsq              :: MVector{16, Float64}             # Buffer to hold <T_i_new T_i_new>
    newState            :: Vector{Complex{Float64}}         # Buffer to hold potential new local state
end

#Initialization from unitcell generated with LatticePhysics
function Configuration( 
    latticename          :: String,                    #Lattice name
    L                    :: Int64,                     #Linear lattice size
    uc                   :: Unitcell{Site{Int64, D}, Bond{Int64, D}} where D, #Unitcell generated with LatticePhysics
    interactions         :: Vector{Matrix{Float64}},       #Vector containing different interactions associated to bond labels in uc
    onsiteInteraction    :: Vector{Float64},      #Onsite interactions
    n                    :: Int64                      #Partons per site (n = 1 -> fundamental or n = 2 -> self-conjugate representation);
    )                    :: Configuration

    #Generate lattice
    l = getLatticePeriodic(uc, L)
    nSites = length(l.sites)

    unitVectors = l.unitcell.lattice_vectors
    basis = point.(uc.sites)
    positions = point.(l.sites)
    basisLabels = label.(l.sites)

    #Read of interactions from bonds
    orgBonds = organizedBondsFrom(l)
    interactionSites = [to.(bond) for bond in orgBonds]
    interactionLabels = [label.(bond) for bond in orgBonds]

    #Generators dependent on filling
    generators = getGenerators(n)
    generatorsSq = generators .^2

    #dimension of local Hilbertspace
    d = length(generators[1][:, 1])

    #Generate random initial state
    state = [getRandomState(d) for _ in 1:length(positions)]

    #Calculate initial spin expectations 
    spinExpectation = [computeSpinExpectation(state[i], generators) for i in 1:nSites]
    spinSqExpectation = [computeSpinExpectation(state[i], generatorsSq) for i in 1:nSites]
    
    #Initialize buffers
    dT = zeros(Float64, length(generators))
    newT = zeros(Float64, length(generators))
    newTsq = zeros(Float64, length(generators))
    newState = zeros(Complex{Float64}, d)
    
    Configuration(latticename, L, nSites, unitVectors, basis, positions, basisLabels, interactionSites,
                  interactionLabels, interactions, onsiteInteraction, d, generators, generatorsSq, 
                  state, spinExpectation, spinSqExpectation, dT, newT, newTsq, newState)
end

##Get functions (using views for slices of arrays)
function Base.length(cfg :: Configuration) :: Int64
    return cfg.nSites
end

@inline function getPosition(cfg :: Configuration, i :: Int) :: Vector{Float64}
    return cfg.positions[i]
end

@inline function getBasis(cfg :: Configuration) :: Vector{Vector{Float64}}
    return cfg.basis
end

@inline function getInteractionSites(cfg :: Configuration, i :: Int) :: Vector{Int64}
    return cfg.interactionSites[i]
end

@inline function getInteractionLabels(cfg :: Configuration, i :: Int) :: Vector{Int64}
    return cfg.interactionLabels[i]
end

@inline function getInteraction(cfg :: Configuration, label :: Int) :: SMatrix{16, 16, Float64, 256}
    return cfg.interactions[label]
end

@inline function getOnsiteInteraction(cfg :: Configuration) :: SVector{16, Float64}
    return cfg.onsiteInteraction
end

@inline function getGenerators(cfg :: Configuration) :: Vector{Matrix{Complex{Float64}}}
    return cfg.generators
end

@inline function getGeneratorsSq(cfg :: Configuration) :: Vector{Matrix{Complex{Float64}}}
    return cfg.generatorsSq
end

@inline function getState(cfg :: Configuration, i :: Int64) :: Vector{ComplexF64}
    return cfg.state[i]
end

@inline function getSpinExpectation(cfg :: Configuration, i :: Int64) :: MVector{16, Float64}
    return cfg.spinExpectation[i]
end

@inline function getSpinSqExpectation(cfg :: Configuration, i :: Int64) :: MVector{16, Float64}
    return cfg.spinSqExpectation[i]
end

@inline function getSpinExpectation(cfg :: Configuration) :: Vector{MVector{16, Float64}} 
    return cfg.spinExpectation
end

@inline function getSpinSqExpectation(cfg :: Configuration) :: Vector{MVector{16, Float64}} 
    return cfg.spinSqExpectation
end

#Compute spin expectation given a state and generators
function computeSpinExpectation(state :: AbstractVector{Complex{Float64}}, generators :: Vector{Matrix{Complex{Float64}}}) :: Vector{Float64}
    return [real(dot(state, generators[mu], state)) for mu in 1:length(generators)]
end

#Compute spin expectation in-place without allocations
function computeSpinExpectation!(newT, state :: AbstractVector{Complex{Float64}}, generators :: Vector{Matrix{Complex{Float64}}}) :: Nothing
    @inbounds for i in eachindex(newT)
        newT[i] = real(dot(state, generators[i], state))
    end
    return nothing
end

#Calculate exchange energy
function exchangeEnergy(
    T_i :: AbstractVector,
    M  :: AbstractMatrix,
    T_j :: AbstractVector
    )  :: Float64

    val = 0.0

    @turbo for k in eachindex(T_i) 
        for l in eachindex(T_j)
            val += T_i[k] * M[k, l] * T_j[l]
        end 
    end 

    return val 
end

#Calculate onsite energy
function onsiteEnergy(
    Tsq :: AbstractVector,
    M   :: AbstractVector
    )   :: Float64

    val = 0.0

    @turbo for k in eachindex(Tsq)
        val += Tsq[k] * M[k]
    end

    return val
end

#Get full energy of lattice
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

#Get energy difference for new spin expectation newT, and newTSq
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
    @turbo cfg.dT .= newT .- T

    for j in eachindex(interactionSites)
        T_j = getSpinExpectation(cfg, interactionSites[j])
        M = getInteraction(cfg, interactionLabels[j])
        dE += exchangeEnergy(cfg.dT, M, T_j)
    end
    return dE
end

#Get energy difference for new spin expectation newT, and newTSq, in-place
function getEnergyDifference!(cfg :: Configuration, i :: Int64) :: Float64
    dE = 0.0
    #Onsite Interaction
    Tsq = getSpinSqExpectation(cfg, i)
    dE += onsiteEnergy(cfg.newTsq, getOnsiteInteraction(cfg)) - 
          onsiteEnergy(Tsq, getOnsiteInteraction(cfg))

    #Exchange Interaction
    T = getSpinExpectation(cfg, i)
    @turbo cfg.dT .= cfg.newT .- T

    interactionSites  = getInteractionSites(cfg, i)
    interactionLabels = getInteractionLabels(cfg, i)

    for j in eachindex(interactionSites)
        T_j = getSpinExpectation(cfg, interactionSites[j])
        M = getInteraction(cfg, interactionLabels[j])
        dE += exchangeEnergy(cfg.dT, M, T_j)
    end
    return dE
end