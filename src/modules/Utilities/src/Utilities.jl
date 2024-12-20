module Utilities
    export solvation_free_energy_with_overlap_penalty
    export get_flat_realization, get_matrix_realization_per_mol, get_point3f_realization, get_point_vector_realization
    export get_matched_distances_between_transformation_offsets, average_offset_distance
    export sum_of_permutation
    export get_initial_state
    export TEMPLATES, TWOTMVSU_EXPERIMENTAL_ASSEMBLY

    using Distances
    using Distributions
    using Rotations
    using StaticArrays
    using GeometryBasics

    include("realization.jl")
    include("configuration_distances.jl")
    include("initialization.jl")
    include("template_data/asymmetric_unit_templates.jl")
    include("template_data/2su_experimental_assembly.jl")
end