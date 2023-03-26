#Interaction Matrix
#3x# matrix for fast products
struct InteractionMatrix
    m11::Float64
    m12::Float64
    m13::Float64
    m21::Float64
    m22::Float64
    m23::Float64
    m31::Float64
    m32::Float64
    m33::Float64
end

#Empty interactionMatrix
function InteractionMatrix()
    return InteractionMatrix(zeros(Float64,9)...)
end

#From normal matrix
function InteractionMatrix(M::AbstractMatrix)
    @assert size(M) == (3,3) "Interaction matrix must be of size 3x3"
    return InteractionMatrix(M[1,1], M[1,2], M[1,3], M[2,1], M[2,2], M[2,3], M[3,1], M[3,2], M[3,3])
end

function InteractionToMatrix(int :: InteractionMatrix) :: Matrix{Float64}
    return [int.m11 int.m12 int.m13;
            int.m21 int.m22 int.m23;
            int.m31 int.m32 int.m33;]
end

#Fast generalized dot product (s_i, M s_j) without allocations
function LinearAlgebra.dot(s_i :: AbstractVector, M :: InteractionMatrix, s_j :: AbstractVector)
    return @fastmath @inbounds s_i[1] * (M.m11 * s_j[1] + M.m12 * s_j[2] + M.m13 * s_j[3]) + s_i[2] * (M.m21 * s_j[1] + M.m22 * s_j[2] + M.m23 * s_j[3]) + s_i[3] * (M.m31 * s_j[1] + M.m32 * s_j[2] + M.m33 * s_j[3])
end

# Diagonal 3x3 fast interaction matrix
struct DiagInteractionMatrix
    m11 :: Float64
    m22 :: Float64
    m33 :: Float64
end

function DiagInteractionMatrix(V :: AbstractVector)
    @assert length(V) == 3 "Interaction matrix must be of size 3x3"
    return DiagInteractionMatrix(V[1], V[2], V[3])
end

function InteractionToVector(int :: DiagInteractionMatrix) :: Vector{Float64}
    return [int.m11, int.m22, int.m33]
end

function LinearAlgebra.dot(s_i :: AbstractVector{Float64}, M :: DiagInteractionMatrix, s_j :: AbstractVector{Float64})
    return @fastmath @inbounds s_i[1] * M.m11 * s_j[1] + s_i[2] * M.m22 * s_j[2] + s_i[3] * M.m33 * s_j[3] 
end

function LinearAlgebra.dot(Tsq :: AbstractVector{Float64}, M :: DiagInteractionMatrix)
    return @fastmath @inbounds Tsq[1] * M.m11 + Tsq[2] * M.m22 + Tsq[3] * M.m33
end

function mult(a :: Float64, M :: DiagInteractionMatrix) :: DiagInteractionMatrix
    return DiagInteractionMatrix(M.m11 * a, M.m22 * a, M.m33 * a)
end
