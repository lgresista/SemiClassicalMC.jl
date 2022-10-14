#Interface function to initialize configuration for different models
function initializeCfg(
        model       :: String,
        latticename :: String,
        J           :: AbstractVector,
        L           :: Int;
        n = 2       :: Int
        )           :: Configuration
        
    if model == "nearest-neighbour"
        initializeCfgNN(J, latticename, L; n = n)
    elseif model == "heisenberg-xxz"
        initializeCfgHeisenbergxxz(J, latticename, L; n = n)
    elseif model == "su2xsu2"
        initializeCfgSU2xSU2(J, latticename, L; n = n)
    elseif model == "tg-hbn"
        initializeCfgTgHbn(J, L; n = n)
    else
        @error "model $model not implemented"
    end
end

## Initialization of specific models
## General NN model for lattice with equivalent nn sites
function initializeCfgNN(
        J           :: Vector{Vector{Matrix{Float64}}},
        latticename :: String,
        L           :: Int;
        n           :: Int = 2
)                   :: Configuration

        @assert length(J) == 2 "Need on-site and neareast neighbour interaction"
        
        #onsite Interaction:
        onsiteInteraction = DiagInteraction(J[1])
        interactions = [Interaction(J[2])]

        # Initialize lattice with bonds (nn bonds are automatically added)
        uc = getUnitcell(Symbol(latticename), Int64, Int64)
        l = getLatticePeriodic(uc, L) 
        orgBonds = organizedBondsFrom(l)

        # Site labels for sublattice magnetization
        siteLabels = getSiteLabels(l, latticename)

        #get interactionSites and interaction labels
        interactionSites = [to.(bond) for bond in orgBonds] #nn bonds
        interactionLabels = [label.(bond) for bond in orgBonds] #all nn bonds have same label = 1

        #Site positions
        positions = point.(l.sites)

        #Generators
        generators = getGenerators(n)

        #Random initial state
        d = size(generators[1], 1)
        state = hcat([getRandomState(d) for _ in 1:length(positions)]...)

        return Configuration(latticename, L, uc.lattice_vectors, point.(uc.sites), positions, siteLabels, interactionSites, interactionLabels, interactions, onsiteInteraction, generators, state)
end

# SU2xSU2
function initializeCfgSU2xSU2(
    J           :: Vector{<:Number},
    latticename :: String,
    L           :: Int;
    n = 2       :: Int64,
    )           :: Configuration

    @assert length(J) == 3 "Need J = [J, Js, Jv]"
    
    J, Js, Jv = J

    onsiteInteraction = getZeroDiagInteraction()
    interactions = [Interaction(
    diagm(ones(3)),

    [J       0.0      0.0
    0.0      J        0.0
    0.0      0.0      J],

    diagm(ones(3)) * Js,
    
    [Jv      0.0      0.0
    0.0      Jv       0.0
    0.0      0.0      Jv])
    ]

    ## Initialize lattice with bonds
    uc = getUnitcell(Symbol(latticename), Int64, Int64)
    l = getLatticePeriodic(uc, L)
    orgBonds = organizedBondsFrom(l)

    siteLabels = getSiteLabels(l, latticename)

    #get interactionSites and interaction labels
    interactionSites = [to.(bond) for bond in orgBonds]
    interactionLabels = [label.(bond) for bond in orgBonds]

    #Site positions
    positions = point.(l.sites)

    #Generate random configuration
    generators = getGenerators(n)

    #Random state
    d = size(generators[1], 1)
    state = hcat([getRandomState(d) for _ in 1:length(positions)]...)
  
    return Configuration(latticename, L, uc.lattice_vectors, point.(uc.sites), positions, siteLabels, interactionSites, interactionLabels, interactions, onsiteInteraction, generators, state)
end


## Heisenberg-XXZ
function initializeCfgHeisenbergxxz(
        J           :: Vector{<:Number},
        latticename :: String,
        L           :: Int;
        n = 2       :: Int64,
        )           :: Configuration

    Jx, Jz, Js, Jvx, Jvz = J

    onsiteInteraction = getZeroDiagInteraction()

    interactions = [Interaction(
    diagm(ones(3)),    

    [Jx      0.0      0.0
    0.0      Jx       0.0
    0.0      0.0      Jz],    

    diagm(ones(3)) * Js,    

    [Jvx    0.0    0.0
    0.0    Jvx    0.0
    0.0    0.0    Jvz]
    )]

    ## Initialize lattice with bonds
    uc = getUnitcell(Symbol(latticename), Int64, Int64)
    l = getLatticePeriodic(uc, L)
    orgBonds = organizedBondsFrom(l)

    siteLabels = getSiteLabels(l, latticename)

    #get interactionSites and interaction labels
    interactionSites = [to.(bond) for bond in orgBonds]
    interactionLabels = [label.(bond) for bond in orgBonds]

    #Site positions
    positions = point.(l.sites)

    #Generate random configuration
    generators = getGenerators(n)
    d = length(generators[1][:, 1])
    state = hcat([getRandomState(d) for _ in 1:length(positions)]...)

    return Configuration(latticename, L, uc.lattice_vectors, point.(uc.sites), positions, siteLabels, interactionSites, interactionLabels, interactions, onsiteInteraction, generators, state)
end


## TG/h-BN (SU2xU1 symmetric)
function initializeCfgTgHbn(
    J     :: Vector{<:Number},
    L     :: Int;
    n = 2 :: Int64,
    )    :: Configuration

    @assert length(J) == 5 "Need J = [J1, J2, Jp1, Jp2, JH]"
    J1, J2, Jp1, Jp2, JH = J/4 #Normalize by four as pfFRG normalizes by 1/8 but has ij and ji bonds

    #Define interaction matrices
    onsiteInteraction = DiagInteraction([1.0, 1.0, 1.0], 
                                        [0.0, 0.0, JH/2], 
                                        [-JH/2, -JH/2, -JH/2],
                                        [0.0, 0.0, JH/2]
                                        )

    nn1Interaction = Interaction(
    diagm(ones(3)) * 1.0,

    [J1 + Jp1 Jp2      0.0
    -Jp2      J1 + Jp1 0.0
    0.0      0.0      J1],

    diagm(ones(3)) * J1,

    [J1 + Jp1 Jp2      0.0
    -Jp2      J1 + Jp1 0.0
    0.0      0.0      J1]
    )

    nn2Interaction = Interaction(
    diagm(ones(3)) * 1.0,

    [J1+Jp1  -Jp2      0.0
    +Jp2      J1+Jp1   0.0
    0.0      0.0      J1],

    diagm(ones(3)) * J1,

    [J1 + Jp1 -Jp2      0.0
    +Jp2      J1 + Jp1  0.0
    0.0      0.0       J1]
    )

    nnnInteraction = Interaction(
    diagm(ones(3)) * 1.0,
    Matrix{Float64}(I, 3, 3) * J2,
    diagm(ones(3)) * J2,
    Matrix{Float64}(I, 3, 3) * J2,
    )

    interactions = [nn1Interaction, nn2Interaction, nnnInteraction]

    ## Initialize lattice with bonds
    uc = getUnitcell(:triangular, Int64, Int64)
    uc.bonds = [] #clear existing bonds

    #Nearest neighbor bonds with labels breaking C6 symmetry
    #Label = 1
    addBond!(uc, 1, 1, 1, (1, 0),  false)
    addBond!(uc, 1, 1, 1, (0, -1), false)
    addBond!(uc, 1, 1, 1, (-1, 1), false)
    #Label = 2
    addBond!(uc, 1, 1, 2, (-1, 0), false)
    addBond!(uc, 1, 1, 2, (0, 1),  false)
    addBond!(uc, 1, 1, 2, (1, -1), false)

    #Symmetric next-nearest neighbor bonds
    #Label = 3
    addBond!(uc, 1, 1, 3, (1, 1),  true)
    addBond!(uc, 1, 1, 3, (-1, 2), true)
    addBond!(uc, 1, 1, 3, (-2, 1), true)

    l = getLatticePeriodic(uc, L)
    orgBonds = organizedBondsFrom(l)
    siteLabels = getSiteLabels(l, "triangular")

    #get interactionSites and interaction labels
    interactionSites = [to.(bond) for bond in orgBonds]
    interactionLabels = [label.(bond) for bond in orgBonds]

    #Site positions
    positions = point.(l.sites)

    #Generate random configuration
    generators = getGenerators(n)
    d = length(generators[1][:, 1])
    state = hcat([getRandomState(d) for _ in 1:length(positions)]...)

    return Configuration("triangular", L, uc.lattice_vectors, point.(uc.sites), positions, siteLabels, interactionSites, interactionLabels, interactions, onsiteInteraction, generators, state)
end