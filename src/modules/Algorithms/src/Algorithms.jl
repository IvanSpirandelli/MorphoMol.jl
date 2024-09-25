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

    # Datastrucutres for both hmc and random walk
    struct MorphometricSimulationInput
        template_mol::Matrix{Float64}
        template_radii::Vector{Float64}
        number_of_molecules::Int
        σ_r::Float64
        σ_t::Float64
        probe_radius::Float64
        packing_fraction::Float64
        prefactors::Vector{Float64}
        overlap_jump::Float64
        overlap_slope::Float64
        T::Float64
        ε::Float64
        L::Int
    end
    
    mutable struct MorphometricSimulationOutput
        states::Vector{Vector{Float64}}
        Es::Vector{Float64}
        Vs::Vector{Float32}
        As::Vector{Float32}
        Ms::Vector{Float32}
        Xs::Vector{Float32}
        OLs::Vector{Float32}
        αs::Vector{Float32}
    end

    mutable struct SimulationOutput
        states::Vector{Vector{Float64}}
        energy_measures::Dict{String, Any}
        algorithm_measures::Dict{String, Any}
    end

    mutable struct SimulationData
        input::MorphometricSimulationInput
        output::MorphometricSimulationOutput
    end

    function add_to_output(x, E, measures, α, output::MorphometricSimulationOutput)
        push!(output.states, deepcopy(x))
        push!(output.Es, E)
        push!(output.Vs, measures[1])
        push!(output.As, measures[2])
        push!(output.Ms, measures[3])
        push!(output.Xs, measures[4])
        push!(output.OLs, measures[5])
        push!(output.αs, α)
    end

    function add_to_output(x, E, energy_measures::Dict{String, Any}, α, output::SimulationOutput)
        push!(output.states, deepcopy(x))
        push!(output.energy_measures["Es"], E)
        for (k, v) in energy_measures
            push!(output.energy_measures[k], v)
        end
        push!(output.algorithm_measures["αs"], α)
    end

    include("hamiltonian_monte_carlo.jl")
    include("random_walk_metropolis.jl")
end