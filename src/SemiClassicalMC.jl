module SemiClassicalMC
        using Random
        using LinearAlgebra
        using LatticePhysics
        using Statistics
        using HDF5
        using Serialization
        using BinningAnalysis
        using LoopVectorization
        using Manopt
        using Manifolds
        using FiniteDifferences
        using ManifoldDiff


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
                initializeCfg,
                getPosition,
                getBasis,
                getInteractionSites,
                getInteractionLabels,
                getInteraction,
                getOnsiteInteraction,
                getGenerators,
                getGeneratorsSq,
                getState,
                getSpinExpectation,
                getSpinSqExpectation,
                getSpinExpectation,
                getSpinSqExpectation,
                computeSpinExpectation,
                computeSpinExpectation!,
                getEnergy,
                getEnergyDifference,
                getEnergyDifference!
        
        include("Models.jl")
        export  initializeCfg
                
        include("Observables/Observables.jl")
        export  Observables,
                ObservablesGeneric,
                ObservablesTgHbn,
                initializeObservables,
                getMagnetization,
                getCorrelations,
                getRs,
                getProject,
                getKsInBox,
                computeStructureFactor,
                computeStructureFactorVec

        include("IO.jl")
        export  checkpoint!,
                readCheckpoint,
                readCheckpointSA,
                getmeans,
                getenergies

        include("MonteCarlo.jl")
        export  run!,
                runAnnealing!,
                getRandomState,
                getRandomState!,
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
                changeLauncherDir!,
                collectAnnealing
end
