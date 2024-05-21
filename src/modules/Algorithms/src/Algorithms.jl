module Algorithms
    export simulate!, standard_leapfrog!
    export HamiltonianMonteCarlo

    using StaticArrays
    using Rotations

    include("hamiltonian_monte_carlo.jl")
end