
## Basic observables for su4 model with non-specific spin-symmetries

#Finite temperature
mutable struct Observables_generic <: Observables
    energy                    :: ErrorPropagator{Float64}
    magnetizationVec          :: LogBinner{Vector{Float64}}
    correlations              :: LogBinner{Matrix{Float64}}
end

function Observables_generic(cfg :: Configuration)
    return Observables_generic(ErrorPropagator(Float64), LogBinner(zeros(Float64, 15)), LogBinner(zeros(Float64, 15, length(getBasis(cfg)) * length(cfg))))
end

function measure!(
    obs         :: Observables_generic,
    cfg         :: Configuration,
    E           :: Number,
    )           :: Nothing
    push!(obs.energy, E/length(cfg), E^2 / length(cfg)^2)
    #Magnetizations 
    push!(obs.magnetizationVec, getMagnetization(cfg))
    push!(obs.correlations, getCorrelations(cfg))
    return nothing
end

#Save MC means after calculation
function saveMeans!(filename :: String, obs :: Observables_generic, cfg :: Configuration,  β :: Float64) :: Nothing
    h5open(filename, "cw") do file
        file["means/energy/mean"]               = means(obs.energy)[1]
        file["means/energy/error"]              = std_errors(obs.energy)[1]
        file["means/magnetizationVec/mean"]  = mean(obs.magnetizationVec)
        file["means/magnetizationVec/error"] = std_error(obs.magnetizationVec)
        file["means/correlations/mean"]  = mean(obs.correlations)
        file["means/correlations/error"] = std_error(obs.correlations)
        c(e) = β * β * (e[2] - e[1] * e[1]) * length(cfg)
        ∇c(e) = [-2.0 * β * β * e[1] * length(cfg), β * β * length(cfg)]
        heat = mean(obs.energy, c)
        dheat = sqrt(abs(var(obs.energy, ∇c, BinningAnalysis._reliable_level(obs.energy))) / obs.energy.count[BinningAnalysis._reliable_level(obs.energy)])
        file["means/heat/mean"] = heat
        file["means/heat/error"] = dheat
    end
    return nothing
end

#Annealing
mutable struct Observables_generic_annealing <: Observables
    sweep                     :: Vector{Int64}
    β                         :: Vector{Float64}
    energy                    :: Vector{Float64}
    magnetizationVec          :: Vector{Vector{Float64}}
    correlations              :: Vector{Matrix{Float64}}
end

function Observables_generic_annealing(cfg :: Configuration) :: Observables_generic_annealing
    Observables_generic_annealing(Float64[], Float64[], Float64[], Vector{Float64}[], Matrix{Float64}[])
end

function measure!(
    obs         :: Observables_generic_annealing,
    cfg         :: Configuration,
    E           :: Number,
    β           :: Number,
    sweep       :: Int
    )           :: Nothing
    push!(obs.β, β)
    push!(obs.sweep, sweep)
    push!(obs.energy, E/length(cfg))
    push!(obs.magnetizationVec, getMagnetization(cfg))
    push!(obs.correlations, getCorrelations(cfg))
    return nothing
end
