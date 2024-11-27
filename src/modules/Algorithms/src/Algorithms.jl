module Algorithms
    export simulate!, standard_leapfrog!
    export RandomWalkMetropolis
    export SimulatedAnnealing
    export MixedEnergyRandomWalkMetropolis
    export HamiltonianMonteCarlo
    export MorphometricSimulationOutput
    export SimulationData
    export add_to_output

    using Dates
    using LinearAlgebra
    using StaticArrays
    using Rotations

    function add_to_output(measures::Dict{String, Any}, output::Dict{String, Vector})
        for (k, v) in measures
            push!(output[k], v)
        end
    end

    include("hamiltonian_monte_carlo.jl")
    include("random_walk_metropolis.jl")
    include("simulated_annealing.jl")
    include("experimental.jl")
end