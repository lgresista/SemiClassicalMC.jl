#Interface function to initialize configuration for different models
function initializeCfg(
        model       :: String,
        latticename :: String,
        J           :: AbstractVector,
        L           :: Int,
        n           :: Int;
        B           :: Vector{Float64} = zeros(Float64, 15)
        )           :: Configuration
        
    if model == "nearest-neighbor"
        initializeCfgNN(J, latticename, L, n; B = B)
    elseif model == "tg-hbn"
        initializeCfgTgHbn(J, L; n = n, B = B)
    else
        @error "model $model not implemented"
    end
end

#Transform indices (μ,ν,κ,η) to (i, j) for σ^μ τ^

## Initialization of specific models
## General NN model for lattice with equivalent nn sites
function initializeCfgNN(
        J           :: Vector{Matrix{Float64}},
        latticename :: String,
        L           :: Int,
        n           :: Int;
        B           :: Vector{Float64} = zeros(Float64, 15),
)                   :: Configuration

        if length(J) == 1
            onsiteInteraction = zeros(Float64, 15)
            interactions = [J[1]]
        elseif length(J) == 2
            #onsite Interaction:
            onsiteInteraction = diag(J[1])
            interactions = [J[2]]
        else
            @error "Need length(J) = 1 (nn interactions) or length(J) = 2 (on-site and nearest-neighbor interactions)."
        end
  
        # Initialize lattice with bonds (nn bonds are automatically added)
        uc = getUnitcell(Symbol(latticename), Int64, Int64)

        return Configuration(latticename, L, uc, interactions, n; onsiteInteraction = onsiteInteraction, B = B)
end

## TG/h-BN as defined in https://doi.org/10.1103/PhysRevB.99.205150
function initializeCfgTgHbn(
    J :: Vector{<:Number},
    L :: Int;
    n :: Int64 = 2,
    B :: Vector{Float64} = zeros(Float64, 15)     
    )    :: Configuration

    @assert length(J) == 5 "Need J = [J1, J2, Jp1, Jp2, JH]"
    J1, J2, Jp1, Jp2, JH = J/8 #Add 1/8 factor

    ## Onsite interaction -J_H/16 * (1+σσ)(1-τᶻτᶻ) (without constant term)
    onsiteInteraction = [0.0, 0.0, JH/2, -JH/2, 0.0, 0.0, JH/2, -JH/2, 0.0, 0.0, JH/2, -JH/2, 0.0, 0.0, JH/2]

    ##Nearest neighbor interaction
    #Spin coupling
    Js_nn = Matrix{Float64}(I, 4, 4)
    
    #Valley coupling 1 (+DM interaction)
    Jv_nn1 =  [ J1      0.0         0.0        0.0
                0.0     J1 + Jp1    Jp2         0.0
                0.0     -Jp2        J1 + Jp1    0.0
                0.0     0.0         0.0         J1 ]

    #Valley coupling 2 (-DM interaction)
    Jv_nn2 =  [ J1      0.0         0.0        0.0
                0.0     J1 + Jp1    -Jp2         0.0
                0.0     Jp2        J1 + Jp1    0.0
                0.0     0.0         0.0         J1 ]

    J_nn1 = kron(Js_nn, Jv_nn1)[2:end, 2:end] #remove density terms
    J_nn2 = kron(Js_nn, Jv_nn2)[2:end, 2:end] #remove density terms

    ##Next-nearest neighbor (diagonal)
    J_nnn = diagm(fill(J2, 15))

    interactions = [J_nn1, J_nn2, J_nnn]

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

    return Configuration("triangular", L, uc, interactions, n; onsiteInteraction = onsiteInteraction, B = B)
end