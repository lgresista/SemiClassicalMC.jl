
## Basic observables for su4 model with non-specific spin-symmetries
mutable struct ObservablesGeneric <: Observables
    energy                    :: ErrorPropagator{Float64}   #energy per site (will also be used to calculate specific heat)
    σ                         :: LogBinner{Float64}         #|Σᵢ⟨σ_i⟩|
    τ                         :: LogBinner{Float64}         #|Σᵢ⟨τ_i⟩|
    στ                        :: LogBinner{Float64}         #|Σᵢ⟨σ_iτ_i⟩|
    magnetizationVec          :: LogBinner{Vector{Float64}} #|Σᵢ⟨T_i⟩| for every generator T_i
    correlations              :: LogBinner{Matrix{Float64}} #Correlations χ(r_i-r_j), with χ[i] = χ(rs[i])
    rs                        :: Vector{Vector{Float64}}    #All possisble r_i -r_j (including different basis sites)
    project                   :: Matrix{Int64}              #k = project[i,j] maps r(i)-r(j) to r(k) - r(1) (for efficient all-to-all correlations)
end

# Initialize Observables from configuration
function ObservablesGeneric(cfg :: Configuration)
    
    rs = getRs(cfg)
    project = getProject(cfg)

    return ObservablesGeneric(
        ErrorPropagator(Float64), 
        LogBinner(Float64),
        LogBinner(Float64),
        LogBinner(Float64),
        LogBinner(zeros(Float64, 15)),
        LogBinner(zeros(Float64, 15, length(cfg.basis) * length(cfg))),
        rs, project
    )
end

# Function run at each measurement
function measure!(
    obs         :: ObservablesGeneric,
    cfg         :: Configuration,
    E           :: Number,
    )           :: Nothing
    push!(obs.energy, E/length(cfg), E^2 / length(cfg)^2)
    #Magnetizations 
    m = getMagnetization(cfg)
    push!(obs.σ, norm(m[1:3]))
    push!(obs.τ, norm(m[4:6]))
    push!(obs.στ, norm(m[7:15]))
    push!(obs.magnetizationVec, m)
    push!(obs.correlations, getCorrelations(cfg, obs.project))
    return nothing
end

#Means saved after calculation
function saveMeans!(filename :: String, obs :: ObservablesGeneric, cfg :: Configuration,  β :: Float64) :: Nothing
    h5open(filename, "cw") do file
        file["means/energy/mean"]            = means(obs.energy)[1]
        file["means/energy/error"]           = std_errors(obs.energy)[1]
        file["means/σ/mean"]  = mean(obs.σ)
        file["means/σ/error"] = std_error(obs.σ)
        file["means/τ/mean"]  = mean(obs.τ)
        file["means/τ/error"] = std_error(obs.τ)
        file["means/στ/mean"]  = mean(obs.στ)
        file["means/στ/error"] = std_error(obs.στ)
        file["means/magnetizationVec/mean"]  = mean(obs.magnetizationVec)
        file["means/magnetizationVec/error"] = std_error(obs.magnetizationVec)
        file["means/correlations/mean"]      = mean(obs.correlations)
        file["means/correlations/error"]     = std_error(obs.correlations)
        c(e)                                 = β * β * (e[2] - e[1] * e[1]) * length(cfg)
        ∇c(e)                                = [-2.0 * β * β * e[1] * length(cfg), β * β * length(cfg)]
        heat                                 = mean(obs.energy, c)
        dheat                                = sqrt(abs(var(obs.energy, ∇c, BinningAnalysis._reliable_level(obs.energy))) / obs.energy.count[BinningAnalysis._reliable_level(obs.energy)])
        file["means/heat/mean"]              = heat
        file["means/heat/error"]             = dheat
        file["rs"] = vectorToMatrix(obs.rs)
    end
    return nothing
end