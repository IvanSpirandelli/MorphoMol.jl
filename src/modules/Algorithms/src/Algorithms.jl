module Algorithms
    export simulate!
    export HamiltonianMonteCarlo

    using StaticArrays
    using Rotations

    include("hamiltonian_monte_carlo.jl")
end