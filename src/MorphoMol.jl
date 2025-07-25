module MorphoMol
using Rotations
# Modules
include("modules/Algorithms/src/Algorithms.jl")
include("modules/Energies/src/Energies.jl")
include("templates/2su_experimental_assembly.jl")
include("templates/asymmetric_unit_templates.jl")
include("templates/get_templates.jl")
include("configuration_distances.jl")

include("simulation_setup/connected_component_calculations.jl")
include("simulation_setup/energies.jl")
include("simulation_setup/initialization.jl")
include("simulation_setup/perturbation.jl")
include("simulation_setup/realizations.jl")
include("simulation_setup/temperature.jl")

# Tests
include("../tests/Tests.jl")

"""
    convert_flat_state_to_tuples(state::Vector{Float64})

A robust helper function to convert a flat `Vector{Float64}` state representation
into the `Vector{Tuple{QuatRotation, Vector{Float64}}}` representation.
"""
function convert_flat_state_to_tuples(state::Vector{Float64})
    @assert mod(length(state), 6) == 0 "Flat state vector length must be a multiple of 6."
    n_mol = length(state) ÷ 6
    state_tuples = Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}(undef, n_mol)
    for i in 1:n_mol
        idx_start = (i - 1) * 6 + 1
        # Assumes Rotations.jl and QuatRotation are available
        R = QuatRotation(exp(Rotations.RotationVecGenerator(state[idx_start:idx_start+2]...)))
        T = state[idx_start+3:idx_start+5]
        state_tuples[i] = (R, T)
    end
    return state_tuples
end

end #module MorphoMol