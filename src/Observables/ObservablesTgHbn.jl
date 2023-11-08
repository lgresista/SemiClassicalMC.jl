## Observables for the tg-hbn hamiltonian

#Finite temperature
mutable struct ObservablesTgHbn <: Observables
    energy                     :: ErrorPropagator{Float64}
    subMsd                     :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_subMsd              :: ErrorPropagator{Float64}
    subMsx                     :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_subMsx              :: ErrorPropagator{Float64}
    subMsz                     :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_subMsz              :: ErrorPropagator{Float64}
    subMdx                     :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_subMdx              :: ErrorPropagator{Float64}
    subMdz                     :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_subMdz              :: ErrorPropagator{Float64}
    msd                        :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_msd                 :: ErrorPropagator{Float64}
    msx                        :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_msx                 :: ErrorPropagator{Float64}
    msz                        :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_msz                 :: ErrorPropagator{Float64}
    mdx                        :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_mdx                 :: ErrorPropagator{Float64}
    mdz                        :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_mdz                 :: ErrorPropagator{Float64}
    spinChirality              :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_spinChirality       :: ErrorPropagator{Float64}
    valleyChirality            :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_valleyChirality     :: ErrorPropagator{Float64}
    spinValleyChirality        :: LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
    binder_spinValleyChirality :: ErrorPropagator{Float64}
    siteLabels                 :: Vector{Int64}
end

# Initializion
function ObservablesTgHbn(cfg :: Configuration)
    siteLabels = getSiteLabelsTriangular120(cfg)
    return ObservablesTgHbn(ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), 
           LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), 
           LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), LogBinner(Float64), ErrorPropagator(Float64), siteLabels)
end

# Measurement
function measure!(
    obs         :: ObservablesTgHbn,
    cfg         :: Configuration,
    E           :: Number,
    )           :: Nothing
    push!(obs.energy, E/length(cfg), E^2 / length(cfg)^2)
    
    κ_σ = getChirality(cfg,[4, 8, 12])
    push!(obs.spinChirality, κ_σ)
    push!(obs.binder_spinChirality, κ_σ.^2, κ_σ.^4)
    
    κ_τ = getChirality(cfg, 1:3)
    push!(obs.valleyChirality, κ_τ)
    push!(obs.binder_valleyChirality, κ_τ.^2, κ_τ.^4)
    
    κ_στ = getChirality(cfg, 5:7) + getChirality(cfg, 9:11) + getChirality(cfg, 13:15)
    push!(obs.spinValleyChirality, κ_στ)
    push!(obs.binder_spinValleyChirality, κ_στ.^2, κ_στ.^4)

    #Magnetizations (sublattice/full) 
    subM = getSublatticeMagnetization(cfg, obs)
    subMsd = mean([norm(subM[[4, 8, 12], n]) for n in 1:3])
    push!(obs.subMsd,subMsd)
    push!(obs.binder_subMsd,subMsd.^2, subMsd.^4)
    subMdx = mean([norm(subM[[1, 2], n]) for n in 1:3])
    push!(obs.subMdx,subMdx)
    push!(obs.binder_subMdx,subMdx.^2, subMdx.^4)
    subMdz = mean([norm(subM[[3], n]) for n in 1:3])
    push!(obs.subMdz,subMdz)
    push!(obs.binder_subMdz,subMdz.^2, subMdz.^4)
    subMsx = mean([norm(subM[[5, 6, 9, 10, 13, 14], n]) for n in 1:3])
    push!(obs.subMsx,subMsx)
    push!(obs.binder_subMsx,subMsx.^2, subMsx.^4)
    subMsz = mean([norm(subM[[7, 11, 15], n]) for n in 1:3])
    push!(obs.subMsz,subMsz)
    push!(obs.binder_subMsz,subMsz.^2, subMsz.^4)

    m = sum(subM, dims = 2)./size(subM, 2)
    msd = mean(norm(m[[4, 8, 12]]))
    push!(obs.msd, msd)
    push!(obs.binder_msd, msd.^2, msd.^4)
    mdx = mean(norm(m[[1, 2]]))
    push!(obs.mdx, mdx)
    push!(obs.binder_mdx, mdx.^2, mdx.^4)
    mdz = mean(norm(m[[3]]))
    push!(obs.mdz, mdz)
    push!(obs.binder_mdz, mdz.^2, mdz.^4)
    msx = mean(norm(m[[5, 6, 9, 10, 13, 14]]))
    push!(obs.msx, msx)
    push!(obs.binder_msx, msx.^2, msx.^4)
    msz = mean(norm(m[[7, 11, 15]]))
    push!(obs.msz, msz)
    push!(obs.binder_msz, msz.^2, msz.^4)
    return nothing
end

#Save MC means after calculation
function saveMeans!(filename :: String, obs :: ObservablesTgHbn, cfg :: Configuration,  β :: Float64) :: Nothing
    h5open(filename, "cw") do file

        b(v) = 1 - v[2]/(3 * v[1]^2)
        ∇b(v) = [2 * v[2] / (3 * v[1]^3), -1/(3 * v[1]^2)]

        file["means/energy/mean"]  = means(obs.energy)[1]
        file["means/energy/error"] = std_errors(obs.energy)[1]

        file["means/subMsd/mean"]         = mean(obs.subMsd)
        file["means/subMsd/error"]        = std_error(obs.subMsd)
        file["means/binder_subMsd/mean"]  = mean(obs.binder_subMsd, b)
        file["means/binder_subMsd/error"] = sqrt(abs(var(obs.binder_subMsd, ∇b, BinningAnalysis._reliable_level(obs.binder_subMsd))) / obs.binder_subMsd.count[BinningAnalysis._reliable_level(obs.binder_subMsd)])

        file["means/subMsx/mean"]         = mean(obs.subMsx)
        file["means/subMsx/error"]        = std_error(obs.subMsx)
        file["means/binder_subMsx/mean"]  = mean(obs.binder_subMsx, b)
        file["means/binder_subMsx/error"] = sqrt(abs(var(obs.binder_subMsx, ∇b, BinningAnalysis._reliable_level(obs.binder_subMsx))) / obs.binder_subMsx.count[BinningAnalysis._reliable_level(obs.binder_subMsx)])

        file["means/subMsz/mean"]         = mean(obs.subMsz)
        file["means/subMsz/error"]        = std_error(obs.subMsz)
        file["means/binder_subMsz/mean"]  = mean(obs.binder_subMsz, b)
        file["means/binder_subMsz/error"] = sqrt(abs(var(obs.binder_subMsz, ∇b, BinningAnalysis._reliable_level(obs.binder_subMsz))) / obs.binder_subMsz.count[BinningAnalysis._reliable_level(obs.binder_subMsz)])

        file["means/subMdx/mean"]         = mean(obs.subMdx)
        file["means/subMdx/error"]        = std_error(obs.subMdx)
        file["means/binder_subMdx/mean"]  = mean(obs.binder_subMdx, b)
        file["means/binder_subMdx/error"] = sqrt(abs(var(obs.binder_subMdx, ∇b, BinningAnalysis._reliable_level(obs.binder_subMdx))) / obs.binder_subMdx.count[BinningAnalysis._reliable_level(obs.binder_subMdx)])

        file["means/subMdz/mean"]         = mean(obs.subMdz)
        file["means/subMdz/error"]        = std_error(obs.subMdz)
        file["means/binder_subMdz/mean"]  = mean(obs.binder_subMdz, b)
        file["means/binder_subMdz/error"] = sqrt(abs(var(obs.binder_subMdz, ∇b, BinningAnalysis._reliable_level(obs.binder_subMdz))) / obs.binder_subMdz.count[BinningAnalysis._reliable_level(obs.binder_subMdz)])

        file["means/msd/mean"]         = mean(obs.msd)
        file["means/msd/error"]        = std_error(obs.msd)
        file["means/binder_msd/mean"]  = mean(obs.binder_msd, b)
        file["means/binder_msd/error"] = sqrt(abs(var(obs.binder_msd, ∇b, BinningAnalysis._reliable_level(obs.binder_msd))) / obs.binder_msd.count[BinningAnalysis._reliable_level(obs.binder_msd)])

        file["means/msx/mean"]         = mean(obs.msx)
        file["means/msx/error"]        = std_error(obs.msx)
        file["means/binder_msx/mean"]  = mean(obs.binder_msx, b)
        file["means/binder_msx/error"] = sqrt(abs(var(obs.binder_msx, ∇b, BinningAnalysis._reliable_level(obs.binder_msx))) / obs.binder_msx.count[BinningAnalysis._reliable_level(obs.binder_msx)])

        file["means/msz/mean"]         = mean(obs.msz)
        file["means/msz/error"]        = std_error(obs.msz)
        file["means/binder_msz/mean"]  = mean(obs.binder_msz, b)
        file["means/binder_msz/error"] = sqrt(abs(var(obs.binder_msz, ∇b, BinningAnalysis._reliable_level(obs.binder_msz))) / obs.binder_msz.count[BinningAnalysis._reliable_level(obs.binder_msz)])

        file["means/mdx/mean"]         = mean(obs.mdx)
        file["means/mdx/error"]        = std_error(obs.mdx)
        file["means/binder_mdx/mean"]  = mean(obs.binder_mdx, b)
        file["means/binder_mdx/error"] = sqrt(abs(var(obs.binder_mdx, ∇b, BinningAnalysis._reliable_level(obs.binder_mdx))) / obs.binder_mdx.count[BinningAnalysis._reliable_level(obs.binder_mdx)])

        file["means/mdz/mean"]         = mean(obs.mdz)
        file["means/mdz/error"]        = std_error(obs.mdz)
        file["means/binder_mdz/mean"]  = mean(obs.binder_mdz, b)
        file["means/binder_mdz/error"] = sqrt(abs(var(obs.binder_mdz, ∇b, BinningAnalysis._reliable_level(obs.binder_mdz))) / obs.binder_mdz.count[BinningAnalysis._reliable_level(obs.binder_mdz)])

        file["means/spinChirality/mean"]         = mean(obs.spinChirality)
        file["means/spinChirality/error"]        = std_error(obs.spinChirality)
        file["means/binder_spinChirality/mean"]  = mean(obs.binder_spinChirality, b)
        file["means/binder_spinChirality/error"] = sqrt(abs(var(obs.binder_spinChirality,∇b, BinningAnalysis._reliable_level(obs.binder_spinChirality))) / obs.binder_spinChirality.count[BinningAnalysis._reliable_level(obs.binder_spinChirality)])

        file["means/valleyChirality/mean"]         = mean(obs.valleyChirality)
        file["means/valleyChirality/error"]        = std_error(obs.valleyChirality)
        file["means/binder_valleyChirality/mean"]  = mean(obs.binder_valleyChirality, b)
        file["means/binder_valleyChirality/error"] = sqrt(abs(var(obs.binder_valleyChirality, ∇b, BinningAnalysis._reliable_level(obs.binder_valleyChirality))) / obs.binder_valleyChirality.count[BinningAnalysis._reliable_level(obs.binder_valleyChirality)])

        file["means/spinValleyChirality/mean"]         = mean(obs.spinValleyChirality)
        file["means/spinValleyChirality/error"]        = mean(obs.spinValleyChirality)
        file["means/binder_spinValleyChirality/mean"]  = mean(obs.binder_spinValleyChirality, b)
        file["means/binder_spinValleyChirality/error"] = sqrt(abs(var(obs.binder_spinValleyChirality, ∇b, BinningAnalysis._reliable_level(obs.binder_spinValleyChirality))) / obs.binder_spinValleyChirality.count[BinningAnalysis._reliable_level(obs.binder_spinValleyChirality)])


        c(e) = β * β * (e[2] - e[1] * e[1]) * length(cfg)
        ∇c(e) = [-2.0 * β * β * e[1] * length(cfg), β * β * length(cfg)]
        heat = mean(obs.energy, c)
        dheat = sqrt(abs(var(obs.energy, ∇c, BinningAnalysis._reliable_level(obs.energy))) / obs.energy.count[BinningAnalysis._reliable_level(obs.energy)])
        file["means/heat/mean"] = heat
        file["means/heat/error"] = dheat
    end
    return nothing
end

### Functions needed for Measurements

# Save sublattice labels for 120 order
function getSiteLabelsTriangular120(cfg :: Configuration)
    lattice_vectors = normalize.(cfg.unitVectors)
    
    basis = [lattice_vectors[1] + lattice_vectors[2], 2*lattice_vectors[2] - lattice_vectors[1]]
    siteLabels = zeros(length(cfg))

    for i in 1:length(cfg)
        v1 = cfg.positions[i]
        v2 = v1 + lattice_vectors[1]

        sol1 = getPrefactors(v1, basis)
        sol2 = getPrefactors(v2, basis)

        if sol1 .% 1 == [0.0, 0.0]
            siteLabels[i] = 1
        elseif sol2 .% 1 == [0.0, 0.0]
            siteLabels[i] = 2
        else
            siteLabels[i] = 3
        end
    end
    return siteLabels
end

#Claculate sublattice magnetization measuring 120 degree order
function getSublatticeMagnetization(cfg :: Configuration, obs :: ObservablesTgHbn) :: Matrix{Float64}
    labels = sort(unique(obs.siteLabels))
    m      = zeros(15, length(labels))

    for (i, label) in enumerate(labels)
        spins = [getSpinExpectation(cfg, j) for j in 1:length(cfg) if obs.siteLabels[j] == label]
        m[:, i] = sum(spins)/length(spins)
    end

    return m
end

#Get vector chirality in z direction averaged over triangles 
function getChirality(cfg :: Configuration, components :: AbstractArray) :: Float64
    κ = 0.0
    for i in 1:length(cfg)
        #Get left and right pointing triangle from site i 
        i2, i3, i4, i5 = getInteractionSites(cfg, i)[[1, 5, 4, 2]]
        
        #Get corresponding spins
        s1 = getSpinExpectation(cfg, i)[components]
        s2 = getSpinExpectation(cfg, i2)[components]
        s3 = getSpinExpectation(cfg, i3)[components]
        s4 = getSpinExpectation(cfg, i4)[components]
        s5 = getSpinExpectation(cfg, i5)[components]

        #compute staggered spin chirality in z direction
        κ += getChirality(s1, s2, s3) - getChirality(s1, s4, s5)
    end
    return κ/length(cfg)/2
end

#Vector chirality in z direction
function getChirality(v1 :: Vector{Float64}, v2 :: Vector{Float64}, v3 :: Vector{Float64}) :: Float64
    return 2/3/sqrt(3) * (v1 × v2 + v2 × v3 + v3 × v1)[3]
end