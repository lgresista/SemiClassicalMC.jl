# Save launcher for finite temperature Monte Carlo
function saveLauncher(
    path             :: String,
    outfile          :: String,
    model            :: String,
    obstype          :: Type{<:Observables},
    latticename      :: String,
    L                :: Int,
    J                :: AbstractVector,
    T                :: Number,
    T_i              :: Number,
    N_therm          :: Int,
    N_measure        :: Int;
    measurement_rate :: Int = 10, 
    checkpoint_rate  :: Int = 100000,
    report_rate      :: Int = 100000,
    normalize        :: Bool = true,
    overwrite        :: Bool = true,
    n                :: Int = 2 #filling
    )                :: Nothing

    open(path, "w") do file
        write(file, """using SemiClassicalMC

                    launch!(
                        "$outfile",            
                        "$model",
                        $obstype,
                        "$latticename",
                        $L,
                        $J,
                        $T,
                        $T_i,
                        $N_therm,
                        $N_measure;
                        measurement_rate = $measurement_rate, 
                        checkpoint_rate = $checkpoint_rate, 
                        report_rate = $report_rate,
                        normalize = $normalize,
                        overwrite = $overwrite,
                        n = $n
                    )
                    """
        )
    end
    return nothing
end

#Launch finite temperature Monte Carlo with one function
function launch!(
    filename         :: String,
    model            :: String,
    obstype          :: Type{<:Observables},
    latticename      :: String,
    L                :: Int,
    J                :: AbstractVector,
    T                :: Number,
    T_i              :: Number,
    N_therm          :: Int,
    N_measure        :: Int;
    measurement_rate :: Int = 10, 
    checkpoint_rate  :: Int = 100000,
    report_rate      :: Int = 100000,
    overwrite        :: Bool = true,
    normalize        :: Bool = true,
    n                :: Int = 2 #filling
    )                :: Nothing


    if normalize == true
        J /= norm(J)
    end

    if overwrite == true
        println("$filename: Overwrite = true, starting from scratch.")
        cfg = initializeCfg(model, latticename, J, L; n = n)
        obs = initializeObservables(obstype, cfg)
        run!(cfg, obs, T, T_i, N_therm, N_measure, filename; measurement_rate = measurement_rate, checkpoint_rate = checkpoint_rate, report_rate = report_rate);
    else
        println("$filename: Overwrite = false, trying to load data and continue calculations")
        cfg, obs, current_sweep, energy, βs, σ= readCheckpoint(filename, obstype)
        if current_sweep >= N_therm + N_measure
            println("$filename: Already reached $(N_therm + N_measure) sweeps, terminating.")
        else
            run!(cfg, obs, T, T_i, N_therm, N_measure, filename; measurement_rate = measurement_rate, checkpoint_rate = checkpoint_rate, report_rate = report_rate,
                 energy = energy, βs = βs, current_sweep = current_sweep, σ = σ)
        end
    end
    return nothing
end

# Save launcher for simulated annealing + minimization
function saveLauncherAnnealing(
    path             :: String,
    outfile          :: String,
    model            :: String,
    obstype          :: Type{<:Observables},
    latticename      :: String,
    L                :: Int,
    J                :: AbstractVector,
    T_i              :: Number;
    T_fac            :: Number = 0.98,
    N_max            :: Int = 10000000,
    N_o              :: Int = 0,
    N_per_T          :: Int = 10000,
    min_acc_per_site :: Int = 100,
    min_accrate      :: Number = 1e-5,
    measurement_rate :: Int = 1000,
    checkpoint_rate  :: Int = 100000,
    σ_min            :: Float64 = 0.05,
    overwrite        :: Bool = true,
    normalize        :: Bool = true,
    seed             :: Int = abs(rand(Random.RandomDevice(),Int)),
    verbose          :: Bool = false,
    n                :: Int64 = 2
    )                :: Nothing
    
    open(path, "w") do file
        write(file, """using SemiClassicalMC
                    
                    launchAnnealing!(
                        "$outfile",            
                        "$model",           
                        $obstype,
                        "$latticename",         
                        $L,            
                        $J,
                        $T_i;
                        T_fac = $T_fac,  
                        N_max = $N_max,             
                        N_o = $N_o,               
                        N_per_T = $N_per_T,           
                        min_acc_per_site = $min_acc_per_site,  
                        min_accrate = $min_accrate,       
                        measurement_rate = $measurement_rate,  
                        checkpoint_rate = $checkpoint_rate,   
                        σ_min = $σ_min,             
                        overwrite = $overwrite,
                        normalize = $normalize,
                        seed = $seed,              
                        verbose = $verbose,           
                        n = $n              
                    )
                    """
        )
    end
    return nothing
end

# Launch simulated annealing + optimization from one function
function launchAnnealing!(
    filename         :: String,
    model            :: String,
    obstype          :: Type{<:Observables},
    latticename      :: String,
    L                :: Int,
    J                :: Vector{<:Number},
    T_i              :: Number;
    T_fac            :: Number = 0.98,
    N_max            :: Int = 10000,
    N_o              :: Int = 0,
    N_per_T          :: Int = 100,
    min_acc_per_site :: Int = 100,
    min_accrate      :: Number  = 0.0,
    measurement_rate :: Int = 10000,
    checkpoint_rate  :: Int = 10000,
    σ_min            :: Float64 = 0.05,
    overwrite        :: Bool    = true,
    normalize        :: Bool    = true,
    seed             :: Int = abs(rand(Random.RandomDevice(),Int)),
    verbose          :: Bool = false,
    n                :: Int = 2 #filling
    )                :: Nothing
    
    if normalize == true
        J /= norm(J)
    end

    if overwrite == true
        println("$filename: Overwrite = true, starting from scratch.")
        cfg = initializeCfg(model, latticename, J, L; n = n)
        obs = initializeObservables(obstype, cfg)
        runAnnealing!(cfg,
                    obs,
                    T_i,
                    filename; 
                    T_fac = T_fac,
                    N_max = N_max,
                    N_o = N_o,
                    N_per_T = N_per_T,
                    min_acc_per_site = min_acc_per_site,
                    min_accrate = min_accrate,
                    measurement_rate = measurement_rate,
                    checkpoint_rate = checkpoint_rate,
                    σ_min = σ_min,
                    seed = seed,
                    verbose = verbose,
                    )
    else
        println("$filename: Overwrite = false, trying to load data and continue calculations")
        cfg, obs, current_sweep, energy, βs, σ = readCheckpoint(filename, obstype)
        runAnnealing!(cfg,
                    obs,
                    T_i,
                    filename; 
                    T_fac = T_fac,
                    N_max = N_max,
                    N_o = N_o,
                    N_per_T = N_per_T,
                    min_acc_per_site = min_acc_per_site,
                    min_accrate = min_accrate,
                    measurement_rate = measurement_rate,
                    checkpoint_rate = checkpoint_rate,
                    σ_min = σ_min,
                    seed = seed,
                    verbose = verbose,
                    βs = βs,
                    energy = energy,
                    current_sweep = current_sweep,
                    σ = σ
                    )
    end
    return nothing
end


function makeJob(path        :: String,
                 dir         :: String,
                 input       :: String,
                 julia       :: String,
                 sbatch_args :: Dict{String, String}
    )           :: Nothing

    # assert that input is a valid Julia script 
    @assert endswith(input, ".jl") "Input must be *.jl file."

    # make local copy to prevent global modification of sbatch_args
    args = copy(sbatch_args)

    # set thread affinity, if not done already
    if haskey(args, "export")
        if occursin("JULIA_EXCLUSIVE", args["export"]) == false
            args["export"] *= ",JULIA_EXCLUSIVE=1"
        end
    else
        push!(args, "export" => "ALL,JULIA_EXCLUSIVE=1")
    end

    # set working directory, if not done already 
    if haskey(args, "chdir")
        @warn "Overwriting working directory passed via SBATCH dict ..."
        args["chdir"] = dir
    else 
        push!(args, "chdir" => dir)
    end

    # set output file, if not done already
    if haskey(args, "output") == false
        output = split(input, ".jl")[1] * ".out"
        push!(args, "output" => output)
    end

    open(path, "w") do file
        # set SLURM parameters
        write(file, "#!/bin/bash \n")

        for arg in keys(args)
            write(file, "#SBATCH --$(arg)=$(args[arg]) \n")
        end

        write(file, "\n")

        # start calculation
        write(file, "$(julia) -O3 -t \$SLURM_CPUS_PER_TASK $(input)")
    end

    return nothing
end


function makeJobSteps(
    path        :: String,
    dir         :: String,
    input       :: String,
    labels      :: AbstractArray,
    julia       :: String,
    sbatch_args :: Dict{String, String}
)               :: Nothing

    # make local copy to prevent global modification of sbatch_args
    args = copy(sbatch_args)

    # set thread affinity, if not done already
    if haskey(args, "export")
        if occursin("JULIA_EXCLUSIVE", args["export"]) == false
            args["export"] *= ",JULIA_EXCLUSIVE=1"
        end
    else
        push!(args, "export" => "ALL,JULIA_EXCLUSIVE=1")
    end

    # set working directory, if not done already 
    if haskey(args, "chdir")
        @warn "Overwriting working directory passed via SBATCH dict ..."
        args["chdir"] = dir
    else 
        push!(args, "chdir" => dir)
    end

    # set output file, if not done already
    if haskey(args, "output") == false
        output = split(input, ".jl")[1] * ".out"
        push!(args, "output" => output)
    end

    open(path, "w") do file
        # set SLURM parameters
        write(file, "#!/bin/bash \n")

        for arg in keys(args)
            write(file, "#SBATCH --$(arg)=$(args[arg])\n")
        end

        write(file, "\n")
        write(file, "$(julia)\n\n")
        # start calculation
        for i in 1:length(labels)
            write(file, "srun -N1 -n1 --exclusive julia -O3 $(input)_$(labels[i]).jl &\n")
        end
        write(file, "wait\n")
        write(file, "echo all jobsteps finished")
    end

    return nothing
end

function changeLauncher!(filename, changes :: Dict{String, String})
    stream = open(filename, "r")
    launcher = read(stream, String)
    close(stream)

    for (key, value) in changes
        if occursin(key, launcher)
            launcher = replace(launcher, key => value)
        else
            println(key, " does not exist in $filename")
        end
    end

    open(filename, "w") do stream
        write(stream, launcher)
    end
    return nothing
end

function changeLauncherDir!(dir, changes)
    filenames = [name for name in readdir(dir) if endswith(name, ".jl")]
    for filename in filenames
        changeLauncher!(filename, changes)
    end
end