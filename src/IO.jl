
## IO using serialization

#Checkpoint for finite T
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

#read checkpoint
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

#Checkpoint for annealing
function checkpoint!(filename, cfg :: Configuration, sweep :: Int64, βs :: Vector{Float64}, energy :: Vector{Float64}, σ :: Float64)
    h5open(filename, "w") do file
        iobuffer = IOBuffer()

        serialize(iobuffer, cfg)
        file["cfg"] = take!(iobuffer)
        file["current_sweep"] = sweep
        file["energy"] = energy
        file["βs"] = βs
        file["σ"] = σ
    end
    return nothing
end

#read checkpoint
function readCheckpointSA(filename) :: Tuple{Configuration, Int64, Vector{Float64}, Vector{Float64}, Float64}
    h5open(filename) do file
        cfg = deserialize(IOBuffer(read(file["cfg"])))
        current_sweep = read(file["current_sweep"])
        energy = read(file["energy"])
        βs = read(file["βs"])
        σ = read(file["σ"])
        return cfg, current_sweep, energy, βs, σ
    end
end

#Collect means for files with diferent temperatures
function getmeans(
    filename,  # function filename(T) = "path-to-file-for-temperature-T"
    outfile,   # filename where means get stored
    Ts         # Temperatures for which to calculate means
)
    #initialize mean dictionary
    println("Initializing output"); flush(stdout)
    means = Dict()
    idxs = Dict()
    f = h5open(filename(Ts[1]))["means"]
    mkeys = keys(f)
    for key in mkeys
        s = size(f[key * "/mean"])
        if length(s) > 1
            t = typeof(read(f[key * "/mean"])).parameters[1]
        else
            t = typeof(read(f[key * "/mean"]))
        end

        means[key * "/mean"] = zeros(t, s..., length(Ts))
        means[key * "/error"] = zeros(t, s..., length(Ts))
        idxs[key] = s
    end
    close(f)

    println("Collecting means for $(length(Ts)) temperatures, this may take a while..."); flush(stdout)
    #collect means
    for i in eachindex(Ts)
        print("$i|"); flush(stdout)
        try
            h5open(filename(Ts[i])) do f
                for key in mkeys
                    key
                    idx = CartesianIndices((idxs[key]..., i:i))
                    means[key * "/mean"][idx] .= read(f["means/" * key * "/mean"])
                    means[key * "/error"][idx] .= read(f["means/" * key * "/error"])
                end
            end
        catch
            println("ERROR: Could not load means from ", filename(Ts[i]));flush(stdout)
        end
    end
    println("\nSaving means"); flush(stdout)
    #save means to file
    h5open(outfile, "w") do f
        f["Ts"] = collect(Ts)
        for key in keys(means)
            f[key] = means[key]
        end
    end
    println("Done")
end

#Get energy measurements (for histogram analysis) for multiple temperatures (usually around T_c)
function getenergies(filename, outfile, Ts)
    println("Collecting energies for $(length(Ts)) temperatures, this may take a while..."); flush(stdout)

    energies = Array{Vector{Float64}}(undef, length(Ts))
    
    #collect means
    for i in eachindex(Ts)
        print("$i|"); flush(stdout)
        try
            h5open(filename(Ts[i])) do f
                energies[i] = read(f["energy"])
            end
        catch
            energies[i] = zeros(1)
            println("ERROR: Could not load measurements from ", filename(Ts[j], labels[i]));flush(stdout)
        end
    end
    println("\nSaving energies"); flush(stdout)
    #save means to file
    h5open(outfile, "w") do f
        for i in eachindex(Ts)
            f["$(Ts[i])"] = energies[i]
        end    
    end
    println("Done")
end