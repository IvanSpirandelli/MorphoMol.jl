module Utilities
    export solvation_free_energy_with_overlap_penalty
    export state_to_poly, poly_to_state
    export get_flat_realization
    export get_matched_distances_between_transformation_offsets, average_offset_distance
    export sum_of_permutation

    using StaticArrays
    using Rotations
    using Distances

    include("poly_read_write.jl")
    include("realization.jl")
    include("configuration_distances.jl")
end