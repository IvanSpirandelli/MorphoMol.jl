module Algorithms
    export simulate!, standard_leapfrog!
    export HamiltonianMonteCarlo
    export MorphometricSimulationOutput
    export SimulationData
    export RandomWalkMetropolis

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
end