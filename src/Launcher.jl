# Save launcher for finite temperature Monte Carlo and predefined models
function saveLauncher(
    filename            :: String,              # Path/to/launcherfile.jl
    model               :: String,              # E.g. "nearest-neighbor" or "tg-hbn" or custom
    n                   :: Int,                 # Filling
    obstype             :: Type{<:Observables}, # Observable subtype determining what is measured (e.g. ObservablesGeneric)
    latticename         :: String,              # E.g. triangular, honeycomb, ... (check LatticePhysics)
    L                   :: Int,                 # Linear lattice size
    J                   :: AbstractVector,      # Couplings given in format required by the specific model
    T                   :: Number,              # Temperature
    T_i                 :: Number,              # Initial temperature for thermalization
    N_therm             :: Int,                 # Number of thermalization sweeps
    N_measure           :: Int;                 # Number of measurement sweeps
    measurement_rate    :: Int = 10,            # Frequency of measurements (in sweeps)
    checkpoint_rate     :: Int = 100000,        # Frequency of checkpoints
    report_rate         :: Int = 100000,        # Frequency of printing progress
    overwrite           :: Bool = true          # Overwrite = false tries to load data from checkpoint
    )                   :: Nothing

    open(filename, "w") do file
        write(file, """using SemiClassicalMC

                    launch!(
                        "$(split(filename, ".jl")[1]).h5",            
                        "$model",
                        $n,
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
                        overwrite = $overwrite
                    )
                    """
        )
    end
    return nothing
end

#Launch finite temperature Monte Carlo with one function
function launch!(
    filename         :: String,              # Path/to/outfile.h5
    model            :: String,              # E.g. "nearest-neighbor" or "tg-hbn" or custom
    n                :: Int,                 #filling 
    obstype          :: Type{<:Observables}, # Observable subtype determining what is measured (e.g. ObservablesGeneric)
    latticename      :: String,              # E.g. triangular, honeycomb, ... (check LatticePhysics)
    L                :: Int,                 # Linear lattice size
    J                :: AbstractVector,      # Couplings given in format required by the specific model
    T                :: Number,              # Temperature at which measurements are taken
    T_i              :: Number,              # Initial temperature for thermalization
    N_therm          :: Int,                 # Number of thermalization sweeps
    N_measure        :: Int;                 # Number of measurement sweeps
    measurement_rate :: Int = 10,            # Frequency of measurements (in sweeps)
    checkpoint_rate  :: Int = 100000,        # Frequency of checkpoints
    report_rate      :: Int = 100000,        # Frequency of printing progress
    overwrite        :: Bool = true,         # Overwrite = false tries to load data from checkpoint
    )                :: Nothing

    if overwrite == true
        #Start from scratch
        println("$filename: Overwrite = true, starting from scratch.")
        cfg = initializeCfg(model, latticename, J, L, n)
        obs = initializeObservables(obstype, cfg)
        run!(cfg, obs, T, T_i, N_therm, N_measure, filename; measurement_rate = measurement_rate, checkpoint_rate = checkpoint_rate, report_rate = report_rate);
    else
        #Try to load checkpoint and continue simulation
        println("$filename: Overwrite = false, trying to load data and continue calculations")
        cfg, obs, current_sweep, energy, βs, σ = readCheckpoint(filename, obstype)
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
    filename         :: String,          # Path/to/outfile.h5
    model            :: String,          # E.g. "nearest-neighbor" or "tg-hbn" or custom
    n                :: Int,             #filling
    latticename      :: String,          # E.g. triangular, honeycomb, ... (check LatticePhysics)
    L                :: Int,             # Linear lattice size
    J                :: AbstractVector,  # Couplings given in format required by the specific model
    T_i              :: Number;          # Initial temperature
    T_fac            :: Number = 0.98,   # Factor by which the temperature gets lowered
    N_max            :: Int = 10000000,  # Maximal number of sweeps before stopping simmulated annealing
    N_o              :: Int = 10,        # Number of optimization sweeps
    N_per_T          :: Int = 10000,     # Maximal number of sweeps before reducing temperature
    min_acc_per_site :: Int = 100,       # If this many updates are accepted per site, reduce temperature
    min_accrate      :: Number = 1e-5,   # Acceptance rate below simulated annealing is stopped
    checkpoint_rate  :: Int = 100000,    # Frequency of checkpoints
    σ_min            :: Float64 = 0.05,  # Minimal stepsize
    seed             :: Int = abs(rand(Random.RandomDevice(),Int)), #seed for RNG
    verbose          :: Bool = false,    #More output in optimization sweeps
    )                :: Nothing
    
    open(filename, "w") do file
        write(file, """using SemiClassicalMC
                    
                    launchAnnealing!(
                        "$(split(filename, ".jl")[1]).h5",            
                        "$model",         
                        $n,  
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
                        checkpoint_rate = $checkpoint_rate,   
                        σ_min = $σ_min,             
                        seed = $seed,              
                        verbose = $verbose
                    )
                    """
        )
    end
    return nothing
end

# Launch simulated annealing + optimization from one function
function launchAnnealing!(
    filename         :: String,         # Path/to/outfile.h5
    model            :: String,         # E.g. "nearest-neighbor" or "tg-hbn" or custom
    n                :: Int,            # filling
    latticename      :: String,         # E.g. triangular,                            honeycomb, ... (check LatticePhysics)
    L                :: Int,            # Linear lattice size
    J                :: AbstractVector, # Couplings given in format required by the specific model
    T_i              :: Number;         # Initial temperature
    T_fac            :: Number = 0.98,  # Factor by which the temperature gets lowered
    N_max            :: Int = 10000,    # Maximal number of sweeps before stopping simmulated annealing
    N_o              :: Int = 0,        # Number of optimization sweeps
    N_per_T          :: Int = 100,      # Maximal number of sweeps before reducing temperature
    min_acc_per_site :: Int = 100,      # If this many updates are accepted per site, reduce temperature
    min_accrate      :: Number  = 0.0,  # Acceptance rate below simulated annealing is stopped
    checkpoint_rate  :: Int = 10000,    # Frequency of checkpoints
    σ_min            :: Float64 = 0.05, # Minimal stepsize
    seed             :: Int = abs(rand(Random.RandomDevice(),Int)), #seed for RNG
    verbose          :: Bool = false    #More output in optimization sweeps
    )                :: Nothing

    cfg = initializeCfg(model, latticename, J, L, n)
    runAnnealing!(cfg,
                T_i,
                filename; 
                T_fac = T_fac,
                N_max = N_max,
                N_o = N_o,
                N_per_T = N_per_T,
                min_acc_per_site = min_acc_per_site,
                min_accrate = min_accrate,
                checkpoint_rate = checkpoint_rate,
                σ_min = σ_min,
                seed = seed,
                verbose = verbose,
                )
 
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
)           :: Nothing


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
            write(file, "srun -N1 -n1 --exclusive julia -O3 -t \$SLURM_CPUS_PER_TASK $(input)_$(labels[i]).jl &\n")
        end
        write(file, "wait\n")
        write(file, "echo all jobsteps finished")
    end

    return nothing
end


function makeJobSteps(path :: String,
    launcherfiles :: Vector{String},
    dir         :: String,
    julia       :: String,
    sbatch_args :: Dict{String, String}
)           :: Nothing


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

    open(path, "w") do file
        # set SLURM parameters
        write(file, "#!/bin/bash \n")

        for arg in keys(args)
            write(file, "#SBATCH --$(arg)=$(args[arg])\n")
        end

        write(file, "\n")
        write(file, "$(julia)\n\n")
        # start calculation
        for fn in launcherfiles
            write(file, "srun -N1 -n1 --exclusive julia -O3 -t \$SLURM_CPUS_PER_TASK $fn &\n")
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