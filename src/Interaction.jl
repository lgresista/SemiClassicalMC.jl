# General spin-valley interaction
struct Interaction
    Js :: InteractionMatrix
    Jv :: InteractionMatrix
    Ks :: InteractionMatrix
    Kv :: InteractionMatrix
end

function Interaction(Js :: Matrix{Float64}, Jv :: Matrix{Float64}, Ks :: Matrix{Float64}, Kv :: Matrix{Float64}) :: Interaction
    return Interaction(InteractionMatrix(Js), InteractionMatrix(Jv), InteractionMatrix(Ks), InteractionMatrix(Kv))
end

function Interaction(J :: Vector{Matrix{Float64}}) :: Interaction
    @assert length(J) == 4 "Interaction needs four matrices: Js, Jv, Ks, Kv"
    return Interaction(J...)
end

function getZeroInteraction() :: Interaction
    return Interaction(zeros(Float64, 3, 3), zeros(Float64, 3, 3), zeros(Float64, 3, 3), zeros(Float64, 3, 3))
end

function unpack(interaction :: Interaction) :: Tuple{InteractionMatrix, InteractionMatrix, InteractionMatrix, InteractionMatrix}
    return interaction.Js, interaction.Jv, interaction.Ks, interaction.Kv
end

function exchangeEnergy(T_i :: Vector{Float64}, M :: Interaction, T_j :: Vector{Float64})
    E = @inbounds @fastmath (dot(T_i[1:3], M.Ks, T_j[1:3])
       + dot(T_i[4:6],   M.Kv, T_j[4:6])
       + dot(T_i[7:9],   M.Jv , T_j[7:9])   * M.Js.m11
       + dot(T_i[10:12], M.Jv , T_j[7:9])   * M.Js.m21
       + dot(T_i[13:15], M.Jv , T_j[7:9])   * M.Js.m31
       + dot(T_i[7:9],   M.Jv , T_j[10:12]) * M.Js.m12
       + dot(T_i[10:12], M.Jv , T_j[10:12]) * M.Js.m22
       + dot(T_i[13:15], M.Jv , T_j[10:12]) * M.Js.m32
       + dot(T_i[7:9],   M.Jv , T_j[13:15]) * M.Js.m13
       + dot(T_i[10:12], M.Jv , T_j[13:15]) * M.Js.m23
       + dot(T_i[13:15], M.Jv , T_j[13:15]) * M.Js.m33
    )
    return E
end

#Diagonal spin-valley interaction (Z^2 x Z^2)
struct DiagInteraction
    Js :: DiagInteractionMatrix
    Jv :: DiagInteractionMatrix
    Ks :: DiagInteractionMatrix
    Kv :: DiagInteractionMatrix
end

function DiagInteraction(Js :: Vector{Float64}, Jv :: Vector{Float64}, Ks :: Vector{Float64}, Kv :: Vector{Float64}) :: DiagInteraction
    return DiagInteraction(DiagInteractionMatrix(Js), DiagInteractionMatrix(Jv), DiagInteractionMatrix(Ks), DiagInteractionMatrix(Kv))
end

function DiagInteraction(J :: Vector{Matrix{Float64}}) :: DiagInteraction
    @assert length(J) == 4 "Interaction needs four matrices: Js, Jv, Ks, Kv"
    for i in eachindex(J)
        @assert size(J[i]) == (3, 3) "Interaction vector needs to contain 3x3 matrices"
        for j in 1:3, k in 1:3
            if j != k && J[i][j, k] != 0.0
                @error "Can not convert offdiagonal matrix to diagonal interaction."
            end
        end
    end
    J_vec = diag.(J)
    return DiagInteraction(J_vec...)
end

function getZeroDiagInteraction() :: DiagInteraction
    return DiagInteraction([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
end

function unpack(interaction :: DiagInteraction) :: Tuple{DiagInteractionMatrix, DiagInteractionMatrix, DiagInteractionMatrix, DiagInteractionMatrix}
    return interaction.Js, interaction.Jv, interaction.Ks, interaction.Kv
end

function onsiteEnergy(Tsq :: Vector{Float64}, M :: DiagInteraction) :: Float64
    E = @inbounds @fastmath (dot(Tsq[1:3], M.Ks)
       + dot(Tsq[4:6], M.Kv)
       + dot(Tsq[7:9] .* M.Js.m11,  M.Jv)
       + dot(Tsq[10:12] .* M.Js.m22,  M.Jv)
       + dot(Tsq[13:15] .* M.Js.m33,  M.Jv)
    )
    return E
end