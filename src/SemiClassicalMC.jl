module SemiClassicalMC
    using Random
    using LinearAlgebra
    using LatticePhysics
    using Statistics
    using HDF5
    using Serialization
    using FiniteDifferences
    using BinningAnalysis
    using Manifolds
    using Manopt
    using PolygonInbounds

  
    include("Generators.jl")
    export  getGenerators,
            expval

    include("InteractionMatrix.jl")

    include("Interaction.jl")
    export  Interaction,
            getZeroInteraction,
            unpack,
            exchangeEnergy,
            DiagInteraction,
            getZeroDiagInteraction,
            onsiteEnergy

    include("Configuration.jl")
    export  Configuration,
            getLatticename,
            getL,
            getUnitVectors,
            getBasis,
            getProject,
            getPosition,
            getSiteLabel,
            getInteractionSites,
            getInteractionLabel,
            getInteraction,
            getOnsiteInteraction,
            getD,
            getGenerators,
            getGeneratorsSq,
            getState,
            getSpinExpectation,
            getSpinSqExpectation,
            computeSpinExpectation,
            getEnergy,
            getEnergyDifference,
            getSiteLabels,
            getSiteLabelsTriangular,
            getPrefactors,
            getProject,
            getRs,
            periodicDistance
    
    include("Models.jl")
    export  initializeCfg

    include("Observables/Observables.jl")
    export  Observables,
            initializeObservables,
            getMagnetization,
            getSublatticeMagnetization,
            getCorrelations,
            getChirality,
            getChirality,
            getCollinearity,
            computeStructureFactor,
            computeStructureFactorVec,
            fold_back!,
            getKsInBox,
            getKsInBZ,
            Observables_tg_hbn_annealing,
            Observables_tg_hbn,
            Observables_generic_annealing,
            Observables_generic,
            measure!,
            saveMeans!

    include("IO.jl")
    export  checkpoint!,
            readCheckpoint,
            saveResult!

    include("MonteCarlo.jl")
    export  runAnnealing!,
            run!,
            getRandomState,
            localUpdate!,
            localOptimization!

    include("Launcher.jl")
    export  saveLauncher,
            launch!,
            saveLauncherAnnealing,
            launchAnnealing!,
            makeJob,
            makeJobSteps,
            changeLauncher!,
            changeLauncherDir!
end
