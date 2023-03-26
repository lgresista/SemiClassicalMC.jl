#Interface function to initialize configuration for different models
function initializeCfg(
        model       :: String,
        latticename :: String,
        J           :: AbstractVector,
        L           :: Int,
        n           :: Int
        )           :: Configuration
        
    if model == "nearest-neighbor"
        initializeCfgNN(J, latticename, L, n)
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
        L           :: Int,
        n           :: Int
)                   :: Configuration

        if length(J) == 1
            onsiteInteraction = getZeroDiagInteraction()
            interactions = [Interaction(J[1])]
        elseif length(J) == 2
            #onsite Interaction:
            onsiteInteraction = DiagInteraction(J[1])
            interactions = [Interaction(J[2])]
        else
            @error "Need length(J) = 1 (nn interactions) or length(J) = 2 (on-site and nearest-neighbor interactions)."
        end
  
        # Initialize lattice with bonds (nn bonds are automatically added)
        uc = getUnitcell(Symbol(latticename), Int64, Int64)

        return Configuration(latticename, L, uc, interactions, onsiteInteraction, n)
end

## TG/h-BN as defined in https://doi.org/10.1103/PhysRevB.99.205150
function initializeCfgTgHbn(
    J     :: Vector{<:Number},
    L     :: Int;
    n = 2 :: Int64,
    )    :: Configuration

    @assert length(J) == 5 "Need J = [J1, J2, Jp1, Jp2, JH]"
    J1, J2, Jp1, Jp2, JH = J/8 #Add 1/8 factor

    #Define interaction matrices
    onsiteInteraction = DiagInteraction([1.0, 1.0, 1.0], 
                                        [0.0, 0.0, JH/2], 
                                        [-JH/2, -JH/2, -JH/2],
                                        [0.0, 0.0, JH/2]
                                        )

    nn1Interaction = Interaction(
    Matrix{Float64}(I, 3, 3) * 1.0,

    [J1 + Jp1 Jp2      0.0
    -Jp2      J1 + Jp1 0.0
    0.0      0.0      J1],

    Matrix{Float64}(I, 3, 3) * J1,

    [J1 + Jp1 Jp2      0.0
    -Jp2      J1 + Jp1 0.0
    0.0      0.0      J1]
    )

    nn2Interaction = Interaction(
    Matrix{Float64}(I, 3, 3) * 1.0,

    [J1+Jp1  -Jp2      0.0
    +Jp2      J1+Jp1   0.0
    0.0      0.0      J1],

    Matrix{Float64}(I, 3, 3) * J1,

    [J1 + Jp1 -Jp2      0.0
    +Jp2      J1 + Jp1  0.0
    0.0      0.0       J1]
    )

    nnnInteraction = Interaction(
    Matrix{Float64}(I, 3, 3) * 1.0,
    Matrix{Float64}(I, 3, 3) * J2,
    Matrix{Float64}(I, 3, 3) * J2,
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

    return Configuration("triangular", L, uc, interactions, onsiteInteraction, n)
end