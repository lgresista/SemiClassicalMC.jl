## IO using serialization
function checkpoint!(filename, cfg :: Configuration, obs :: Observables, sweep :: Int64, βs :: Vector{Float64}, energy :: Vector{Float64}, σ :: Float64)
    h5open(filename, "w") do file
        iobuffer = IOBuffer()

        serialize(iobuffer, cfg)
        file["cfg"] = take!(iobuffer)

        serialize(iobuffer, obs)
        file["obs"] = take!(iobuffer)
        
        file["current_sweep"] = sweep
        file["energy"] = energy
        file["βs"] = βs
        file["σ"] = σ
    end
    return nothing
end

function readCheckpoint(filename, :: Type{Obs}) :: Tuple{Configuration, Obs, Int64, Vector{Float64}, Vector{Float64}, Float64} where {Obs <: Observables}
    h5open(filename) do file
        cfg = deserialize(IOBuffer(read(file["cfg"])))
        obs = deserialize(IOBuffer(read(file["obs"])))
        current_sweep = read(file["current_sweep"])
        energy = read(file["energy"])
        βs = read(file["βs"])
        σ = read(file["σ"])
        return cfg, obs, current_sweep, energy, βs, σ
    end
end

## Save result of minimization (last entry of observables)
function saveResult!( 
    filename :: String,
    obs      :: Observables,
    )        :: Nothing
    h5open(filename, "cw") do f
        res = create_group(f, "res")
        for field in fieldnames(typeof(obs))
            res[String(field)] = getfield(obs, field)[end]
        end
    end
end


############### Helpfer functions to convert between vectors and arrays ###############

function vectorToMatrix(vector :: Vector{Vector{T}}) where T
    ncol = length(vector[1])
    nrow = length(vector)
    return [vector[j][i] for i in 1:ncol, j in 1:nrow]
end

function vectorToArray(vector :: Vector{Matrix{T}}) where T
    n1, n2 = size(vector[1])
    n3 = length(vector)
    return [vector[k][i, j] for i in 1:n1, j in 1:n2, k in 1:n3]
end

function arrayToMatrixVector(array :: Array{T, 3}) where T
    return [array[:, :, k] for k in 1:size(array, 3)]
end

function vectorToArray(vector :: Vector{Vector{Vector{T}}}) where T
    n1 = length(vector[1][1])
    n2 = length(vector[1])
    n3 = length(vector)
    
    return [vector[j][i][h] for h in 1:n1, i in 1:n2, j in 1:n3]
end

function matrixToVector(matrix)
    return [matrix[:, i] for i in 1:size(matrix, 2)]
end

function arrayToVector(array :: Array{T, 3}) where T
    return [[array[:, i, j] for i in 1:size(array, 2)] for j in 1:size(array, 3)]
end