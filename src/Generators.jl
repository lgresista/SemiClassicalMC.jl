#Get generators for a filling of n partons
function getGenerators(n :: Int64) :: Vector{Matrix{Complex{Float64}}}
    #fermionic basis
    basis = [digits(n; base = 2, pad = 4) for n in 0:15]

    #Return S(μ, ν, n, basis) = σ^μ ⊗ τ^ν at filling n
    return [S(0,  0, n, basis), #1
            S(0,  1, n, basis), #2
            S(0,  2, n, basis), #3
            S(0,  3, n, basis), #4
            S(1,  0, n, basis), #5
            S(1,  1, n, basis), #6
            S(1,  2, n, basis), #7
            S(1,  3, n, basis), #8
            S(2,  0, n, basis), #9
            S(2,  1, n, basis), #10
            S(2,  2, n, basis), #11
            S(2,  3, n, basis), #12
            S(3,  0, n, basis), #13
            S(3,  1, n, basis), #14
            S(3,  2, n, basis), #15
            S(3,  3, n, basis)] #16
end

#Function (μ, ν) to i so that getGenerators[get_su4index(μ,ν)] = σ^μ ⊗ σ^ν
function get_su4index(μ :: Int, ν :: Int) :: Int
    return 4*μ + ν + 1
end

#Choosen bassis: Occupation of 
# (up, +), (up, -), (down, +), (down, -)

#spin (s = 1 (up), 2 (down)) and valley(v = 1 (+), 2 (-)) to one index (j = (1, 2, 3, 4) corresponding to the basis above)
function sl_to_j(s, l)
    2s + l - 2
end

#Fermionic operator in matrix representation <a_i|(f_sl = f_j)|a_j> with |a> in basis
function f(s, l, basis)
    mat = zeros(length(basis), length(basis))
    for j in 1:length(basis)
        f_m = copy(basis[j]) # f_{sl}|a_j>
        f_m[sl_to_j(s, l)] -= 1
        
        i = findfirst(basis.==[f_m]) #Use <a_i|a_j> = δ_ij
        if isnothing(i) == false
            mat[i, j] = 1
        end
    end
    return mat
end

#Spin-operator in matrix represenation (for full Hilbert space)
function S_full(μ, ν, basis)
    return sum(f(s, l, basis)'*f(sp, lp, basis) * σ(μ)[s, sp] * σ(ν)[l, lp] for s in 1:2, sp in 1:2, l in 1:2, lp in 1:2)
end

#Projector to physical Hilbert space
function P(O, n, basis)
    projIndices = [sum(b) == n for b in basis]
    return O[projIndices, projIndices]
end

#Spin-operator projected to physical Hilbertspace with filling n
S(μ, ν, n, basis) = P(S_full(μ, ν, basis), n, basis)

#Pauli matrices (including identity matrix) as 3D Array
function get_pauli_matrices(normalization :: Real)
    pauli = zeros(Complex{Float64}, 2, 2, 4)
    pauli[:, :, 1] = [1 0; 0 1]
    pauli[:, :, 2] = [0 1; 1 0]
    pauli[:, :, 3] = [0 -im; im 0]
    pauli[:, :, 4] = [1 0; 0 -1]
    pauli[:, :, 2:4] /= normalization;
    return pauli
end

function σ(μ :: Int; normalization = 1.0) :: Matrix{Complex{Float64}}
    return get_pauli_matrices(normalization)[:, :, μ+1]
end

function expval(M :: AbstractMatrix, state :: AbstractVector) :: Complex{Float64}
    return dot(state, M, state)
end
