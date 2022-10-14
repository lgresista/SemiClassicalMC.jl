
#Finite temperature
mutable struct Observables_tg_hbn <: Observables
    energy                    :: ErrorPropagator{Float64}
    subMsd                    :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    subMsx                    :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    subMsz                    :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    subMdx                    :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    subMdz                    :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    msd                       :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    msx                       :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    msz                       :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    mdx                       :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    mdz                       :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    spinChirality             :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    valleyChirality           :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    spinValleyChirality       :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    correlations              :: LogBinner{Matrix{Float64}}
end

function Observables_tg_hbn(cfg :: Configuration)
    return Observables_tg_hbn(ErrorPropagator(Float64), LogBinner(Float64), LogBinner(Float64), LogBinner(Float64), LogBinner(Float64), 
           LogBinner(Float64), LogBinner(Float64), LogBinner(Float64), LogBinner(Float64), LogBinner(Float64), LogBinner(Float64), 
           LogBinner(Float64), LogBinner(Float64), LogBinner(Float64), LogBinner(zeros(Float64, 15, length(getBasis(cfg)) * length(cfg))))
end

function measure!(
    obs         :: Observables_tg_hbn,
    cfg         :: Configuration,
    E           :: Number,
    )           :: Nothing
    push!(obs.energy, E/length(cfg), E^2 / length(cfg)^2)
    push!(obs.spinChirality, getChirality(cfg, 1:3))
    push!(obs.valleyChirality, getChirality(cfg, 4:6))
    push!(obs.spinValleyChirality, getChirality(cfg, 7:9) + getChirality(cfg, 10:12) + getChirality(cfg, 13:15))

    #Magnetizations (sublattice/full) 
    subM = getSublatticeMagnetization(cfg)
    push!(obs.subMsd, mean([norm(subM[[1, 2, 3], n]) for n in 1:3]))
    push!(obs.subMdx, mean([norm(subM[[4, 5], n]) for n in 1:3]))
    push!(obs.subMdz, mean([norm(subM[[6], n]) for n in 1:3]))
    push!(obs.subMsx, mean([norm(subM[[7, 8, 10, 11, 13, 14], n]) for n in 1:3]))
    push!(obs.subMsz, mean([norm(subM[[9, 12, 15], n]) for n in 1:3]))

    m = sum(subM, dims = 2)./size(subM, 2)
    push!(obs.msd, mean(norm(m[[1, 2, 3]])))
    push!(obs.mdx, mean(norm(m[[4, 5]])))
    push!(obs.mdz, mean(norm(m[[6]])))
    push!(obs.msx, mean(norm(m[[7, 8, 10, 11, 13, 14]])))
    push!(obs.msz, mean(norm(m[[9, 12, 15]])))
    return nothing
end

#Save MC means after calculation
function saveMeans!(filename :: String, obs :: Observables_tg_hbn, cfg :: Configuration,  β :: Float64) :: Nothing
    h5open(filename, "cw") do file
        file["means/energy/mean"]               = means(obs.energy)[1]
        file["means/energy/error"]              = std_errors(obs.energy)[1]
        file["means/subMsd/mean"]               = mean(obs.subMsd)
        file["means/subMsd/error"]              = std_error(obs.subMsd)
        file["means/subMsx/mean"]               = mean(obs.subMsx)
        file["means/subMsx/error"]              = std_error(obs.subMsx)
        file["means/subMsz/mean"]               = mean(obs.subMsz)
        file["means/subMsz/error"]              = std_error(obs.subMsz)
        file["means/subMdx/mean"]               = mean(obs.subMdx)
        file["means/subMdx/error"]              = std_error(obs.subMdx)
        file["means/subMdz/mean"]               = mean(obs.subMdz)
        file["means/subMdz/error"]              = std_error(obs.subMdz)
        file["means/msd/mean"]                  = mean(obs.msd)
        file["means/msd/error"]                 = std_error(obs.msd)
        file["means/msx/mean"]                  = mean(obs.msx)
        file["means/msx/error"]                 = std_error(obs.msx)
        file["means/msz/mean"]                  = mean(obs.msz)
        file["means/msz/error"]                 = std_error(obs.msz)
        file["means/mdx/mean"]                  = mean(obs.mdx)
        file["means/mdx/error"]                 = std_error(obs.mdx)
        file["means/mdz/mean"]                  = mean(obs.mdz)
        file["means/mdz/error"]                 = std_error(obs.mdz)
        file["means/spinChirality/mean"]        = mean(obs.spinChirality)
        file["means/spinChirality/error"]       = std_error(obs.spinChirality)
        file["means/valleyChirality/mean"]      = mean(obs.valleyChirality)
        file["means/valleyChirality/error"]     = std_error(obs.valleyChirality)
        file["means/spinValleyChirality/mean"]  = mean(obs.spinValleyChirality)
        file["means/spinValleyChirality/error"] = mean(obs.spinValleyChirality)
        file["means/correlations/mean"]         = mean(obs.correlations)         
        file["means/correlations/error"]        = std_error(obs.correlations)

        c(e) = β * β * (e[2] - e[1] * e[1]) * length(cfg)
        ∇c(e) = [-2.0 * β * β * e[1] * length(cfg), β * β * length(cfg)]
        heat = mean(obs.energy, c)
        dheat = sqrt(abs(var(obs.energy, ∇c, BinningAnalysis._reliable_level(obs.energy))) / obs.energy.count[BinningAnalysis._reliable_level(obs.energy)])
        file["means/heat/mean"] = heat
        file["means/heat/error"] = dheat
    end
    return nothing
end

## Observables for the tg-hbn hamiltonian
#Annealing
mutable struct Observables_tg_hbn_annealing <: Observables
    sweep                     :: Vector{Int64}
    β                         :: Vector{Float64}
    energy                    :: Vector{Float64}
    subMagnetizationVec       :: Vector{Matrix{Float64}}
    correlations              :: Vector{Matrix{Float64}}
    spinChirality             :: Vector{Float64}
    valleyChirality           :: Vector{Float64}
    spinValleyChirality       :: Vector{Float64}
end

function Observables_tg_hbn_annealing(cfg :: Configuration) :: Observables_tg_hbn_annealing
    Observables_tg_hbn_annealing(Float64[], Float64[], Float64[], Matrix{Float64}[], Matrix{Float64}[],
                    Float64[], Float64[], Float64[])
end

function measure!(
    obs         :: Observables_tg_hbn_annealing,
    cfg         :: Configuration,
    E           :: Number,
    β           :: Number,
    sweep       :: Int
    )           :: Nothing
    push!(obs.β, β)
    push!(obs.sweep, sweep)
    push!(obs.energy, E/length(cfg))
    push!(obs.subMagnetizationVec, getSublatticeMagnetization(cfg))
    push!(obs.correlations, getCorrelations(cfg))
    push!(obs.spinChirality, getChirality(cfg, 1:3))
    push!(obs.valleyChirality, getChirality(cfg, 4:6))
    push!(obs.spinValleyChirality, getChirality(cfg, 7:9) + getChirality(cfg, 10:12) + getChirality(cfg, 13:15))
    return nothing
end