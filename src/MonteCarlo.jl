## Different Monte Carlo techniques applicable to general configurations
# Fixed T (with simulated annealing for thermalization)
function run!(
    cfg              :: Configuration,
    obs              :: Observables,
    T                :: Number,
    T_i              :: Number,
    N_therm          :: Int,
    N_measure        :: Int,
    filename         :: String;
    measurement_rate :: Int = 10,
    checkpoint_rate  :: Int = 100000,
    report_rate      :: Int = 100,
    current_sweep    :: Int = 0,
    energy           :: Vector{Float64} = Float64[],
    βs               :: Vector{Float64} = Float64[],
    σ                :: Float64 = 60.0,
    seed             :: Int = abs(rand(Random.RandomDevice(),Int)),
    )                :: Nothing

    #Reseeding RNG
    Random.seed!(seed)
    
    #Initializing output
    println("$filename: Initializaing output ----"); flush(stdout)

    #Initialize betas for thermalization
    N_anneal = ceil(Int64, 3/4 * N_therm)
    Ts = Float64[T_i * (T/T_i)^((n-1)/(N_anneal-1)) for n in 1:N_anneal]
    Ts = [Ts; fill(T, N_therm - N_anneal)]
    βs_therm = 1 ./ Ts

    E = getEnergy(cfg)
    #Thermalization
    if current_sweep <= N_therm 
        println("$filename: Starting thermalization sweeps"); flush(stdout)
    end

    attempted_updates = 0
    accepted_updates = 0
    R = 0.5
    for sweep in current_sweep + 1 : N_therm  ## Go through maximal number of MonteCarlo sweeps
        β = βs_therm[sweep]
        for update in 1 : 2 * length(cfg)
            i = rand(1:length(cfg))
            E, accepted_updates = localUpdate!(cfg, E, accepted_updates, β, i, σ)
            attempted_updates += 1
        end
        
        #Adapt cone width of update
        if sweep % measurement_rate == 0
            push!(βs, β)
            push!(energy, E/length(cfg))
            R = accepted_updates/attempted_updates
            σ = min(σ * 0.5 / (1-R), 60)
            attempted_updates = 0
            accepted_updates = 0
        end 

        if sweep % checkpoint_rate == 0
            checkpoint!(filename, cfg, obs, sweep, βs, energy, σ)
            saveMeans!(filename, obs, cfg, β)
        end

        if sweep % report_rate == 0
            T_round = round(1/β, sigdigits = 4)
            println("$filename: $sweep/$N_therm thermalization sweeps. T = $(T_round), R = $(round(R, sigdigits = 4)), σ = $(round(σ, sigdigits = 4)).")        
        end

    end

    #Measurements
    attempted_updates = 0
    accepted_updates = 0
    β = 1/T

    println("$filename: Starting measurement sweeps at sweep $(max(current_sweep, N_therm)-N_therm)/$N_measure, T = $(round(1/β, sigdigits = 4)). R = $(round(R, sigdigits = 4)), σ = $(round(σ, sigdigits = 4))."); flush(stdout)

    for sweep in max(N_therm + 1, current_sweep + 1) : N_measure + N_therm
        for update in 1 : 2 * length(cfg)
            i = rand(1:length(cfg))
            E, accepted_updates = localUpdate!(cfg, E, accepted_updates, β, i, σ)
            attempted_updates += 1
        end

        if sweep % measurement_rate == 0
            measure!(obs, cfg, E)
            push!(energy, E/length(cfg))
            push!(βs, β)
        end 

        if sweep % checkpoint_rate == 0
            checkpoint!(filename, cfg, obs, sweep, βs, energy, σ)
            saveMeans!(filename, obs, cfg, β)
        end

        if sweep % report_rate == 0
            R = accepted_updates/attempted_updates
            attempted_updates = 0
            accepted_updates = 0
            println("$filename: $(sweep-N_therm)/$N_measure measurement sweeps at $(round(R, sigdigits = 4)) acceptance rate.")
        end
    end

    println("$filename: Measurement sweeps finished. Ending calculation."); flush(stdout)
    
    checkpoint!(filename, cfg, obs, max(current_sweep, N_therm + N_measure), βs, energy, σ)
    saveMeans!(filename, obs, cfg, β)

    println("$filename: Calculation finished"); flush(stdout)
    return nothing
end


# Simulated annealing
function runAnnealing!(
    cfg              :: Configuration,
    T_i              :: Number,
    filename         :: String;
    T_fac            :: Number = 0.98,
    N_max            :: Int = 10000000,
    N_o              :: Int = 0,
    N_per_T          :: Int = 10000,
    min_acc_per_site :: Int = 100,
    min_accrate      :: Number = 1e-5,
    checkpoint_rate  :: Int = 100000,
    σ_min            :: Float64 = 0.05,
    seed             :: Int = abs(rand(Random.RandomDevice(),Int)),
    verbose          :: Bool = false,
    βs               :: Vector{Float64} = Float64[],
    energy           :: Vector{Float64} = Float64[],
    current_sweep    :: Int64 = 1,
    σ                :: Float64 = 60.0
    )                :: Nothing
    
    #Reseed RNG 
    Random.seed!(seed)

    E = getEnergy(cfg)
    
    #Initialize beta
    β = isempty(βs) ? 1/T_i : βs[end]
    β_fac = 1/T_fac
    min_accepted = min_acc_per_site * length(cfg)
    
    R = 1.0     
    sweep = copy(current_sweep)

    if sweep < N_max 
        println("$filename: Starting annealing sweeps at $sweep/$N_max, T = $(round((1/β), sigdigits = 4)), E = $(round(E/length(cfg), sigdigits = 4))"); flush(stdout)
    end

    while sweep < N_max

        attempted_updates = 0
        accepted_updates = 0
        sweep_at_T = 0

        while accepted_updates < min_accepted && sweep_at_T < N_per_T && sweep < N_max
            for update in 1 : 2 * length(cfg)
                idx = rand(1:length(cfg))
                E, accepted_updates = localUpdate!(cfg, E, accepted_updates, β, idx, σ)
                attempted_updates += 1
            end

            push!(energy, E/length(cfg))
            push!(βs, β)

            if sweep % checkpoint_rate == 0
                R_cp = round(accepted_updates/attempted_updates, sigdigits = 4)
                println("$filename: Generating checkpoint after $sweep/$N_max sweeps, \
                T = $(round(1/β, sigdigits = 4)), \
                E = $(round(E/length(cfg), sigdigits = 4)), \
                R = $(round(R_cp, sigdigits = 4)), \
                σ = $(round(σ, sigdigits = 4))"); flush(stdout)
                checkpoint!(filename, cfg, sweep, βs, energy, σ)
            end
            sweep_at_T +=1
            sweep +=1
        end

        #Adapt cone width of update
        R = accepted_updates/attempted_updates
        
        #Check if minimal acceptance rate is reached
        if R > min_accrate
            if verbose
                println("$filename: Reducing temperature at $sweep/$N_max sweeps
                                T = $(round(1/β, sigdigits = 4)),
                                E = $(round(E/length(cfg), sigdigits = 4)),
                                R = $(round(R, sigdigits = 4)),
                                σ = $(round(σ, sigdigits = 4)),
                                accepted updates per site = $(round(accepted_updates/length(cfg), digits = 2)),
                                sweeps = $sweep_at_T")
                        ; flush(stdout)
            end
            #Lower temperature
            β *= β_fac
            #Adapt cone width of update
            σ = max(min(σ * 0.5 / (1-R), 60), σ_min)
        else 
            #End annealing sweeps
            println("$filename: Minimal acceptance rate reached at $sweep/$N_max sweeps, \
                    T = $(round(1/β, sigdigits = 4)), \
                    E = $(round(E/length(cfg), sigdigits = 4)), \
                    R = $(round(R, sigdigits = 4)), \
                    σ = $(round(σ, sigdigits = 4)),
                    accepted updates per site = $(round(accepted_updates/length(cfg), digits = 2)), \
                    sweeps = $sweep_at_T");flush(stdout)   
            break
        end
    end
    
    if sweep == N_max
        println("$filename: Maximum number of annealing sweeps reached at T = $(1/β) ($sweep sweeps, E = $(round(E/length(cfg), sigdigits = 6)))."); flush(stdout)
    end
    
    if N_o > 0
        println("$filename: Starting optimization sweeps"); flush(stdout)
    end

    for o in 1:N_o
        for update in 1 : 2 * length(cfg)
            i = rand(1:length(cfg))
            E = localOptimization!(cfg, E, i)
        end

        push!(energy, E/length(cfg))
        push!(βs, β)

        println("$filename: $o/$N_o optimization sweeps done. E = $(E/length(cfg))"); flush(stdout)
        sweep +=1
    end

    println("$filename: Finished optimization sweeps")
    println("$filename: Saving observables"); flush(stdout)

    checkpoint!(filename, cfg, sweep, βs, energy, σ)

    h5open(filename, "cw") do f
        f["state"] = vectorToMatrix(cfg.state)
    end

    println("$filename: Finished calculation."); flush(stdout)
    return nothing
end



## updates

# Normalized complex vector of dimension d sampled uniformly 
function getRandomState(d :: Int64) :: Vector{Complex{Float64}}
    return normalize(randn(Complex{Float64}, d))
end

function getRandomState(state :: AbstractVector{Complex{Float64}}, σ :: Float64) :: Vector{Complex{Float64}}
    return normalize(state .+ σ .* randn(Complex{Float64}, length(state)))
end

function getRandomState!(newState :: AbstractVector{Complex{Float64}}, state :: AbstractVector{Complex{Float64}}, σ :: Float64) :: Nothing
    newState .= (state .+ σ .* randn(Complex{Float64}, length(state)))
    normalize!(newState)
    return nothing
end

@inline function getRandomState!(
    newState :: StructVector{ComplexF64, NamedTuple{(:re, :im), Tuple{Vector{Float64}, Vector{Float64}}}, Int64}, 
    state    :: StructVector{ComplexF64, NamedTuple{(:re, :im), Tuple{Vector{Float64}, Vector{Float64}}}, Int64}, 
    σ        :: Float64
    )        :: Nothing

    randn!(local_rng(), newState.re); @turbo newState.re .*= σ; @turbo newState.re .+= state.re
    randn!(local_rng(), newState.im); @turbo newState.im .*= σ; @turbo newState.im .+= state.im
    normalize!(newState)
    
    return nothing
end

function localUpdate!(cfg :: Configuration, E :: Float64, accepted_updates :: Int64, β :: Float64, i :: Int64, σ :: Float64) :: Tuple{Float64, Int64}
    
    #Generate new state (stored in cfg.newState)
    state = getState(cfg, i)
    getRandomState!(cfg.newState, state, σ)

    #Calculate <T_i> and <T_i^2> and store in cfg.newT and cfg.newTsq
    computeSpinExpectation!(cfg.newT, cfg.newState, getGenerators(cfg))
    computeSpinExpectation!(cfg.newTsq, cfg.newState, getGeneratorsSq(cfg))

    #Calculate energy difference
    dE = getEnergyDifference!(cfg, i)
    
    #Metropolis update
    p = exp(-β*dE)
    if rand() < p
        cfg.state[i] .= cfg.newState
        @turbo cfg.spinExpectation[i] .= cfg.newT
        @turbo cfg.spinSqExpectation[i] .= cfg.newTsq
        E += dE
        accepted_updates += 1
    end
    return E, accepted_updates
end


function localOptimization!(cfg :: Configuration, E :: Float64, i :: Int64) :: Float64

    M = Sphere(2*cfg.d-1)

    function F(M, realstate)
        #newState = StructArray(realToComplex(realstate))
        newState = StructVector{ComplexF64}(re = realstate[1:cfg.d], im = realstate[cfg.d+1:end])
        newT = computeSpinExpectation(newState, getGenerators(cfg))
        newTsq = computeSpinExpectation(newState, getGeneratorsSq(cfg))
        return getEnergyDifference(cfg, i, newT, newTsq)
    end

    F(realstate) = F(M, realstate)

    r_backend = ManifoldDiff.TangentDiffBackend(ManifoldDiff.FiniteDifferencesBackend()) 

    gradF(M, state) = Manifolds.gradient(M, F, state, r_backend)

    res = gradient_descent(
            M,
            F,
            gradF,
            [getState(cfg, i).re; getState(cfg, i).im]
            ;
            stepsize=ArmijoLinesearch(initial_stepsize = 1.0, retraction_method = ExponentialRetraction(), contraction_factor = 0.99, sufficient_decrease = 0.5),
            stopping_criterion = (StopWhenAny(StopAfterIteration(200),StopWhenGradientNormLess(1e-3))),
            #debug=[
            #    :Iteration,
            #    (:Change, "change: %1.9f | "),
            #    (:Cost, " F(x): %1.11f | "),
            #    "\n",
            #    :Stop,
            #],
            )
    E += F(S, res)
    newState = StructVector{ComplexF64}(re = res[1:cfg.d], im = res[cfg.d+1:end])
    cfg.state[i] .= newState
    cfg.spinExpectation[i] .= computeSpinExpectation(newState, getGenerators(cfg))
    cfg.spinSqExpectation[i] .= computeSpinExpectation(newState, getGeneratorsSq(cfg))
    return E
end