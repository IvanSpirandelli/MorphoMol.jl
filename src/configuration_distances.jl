using Combinatorics
using LinearAlgebra
using Statistics

"""
    get_configuration_distance(state_a, state_b, templates_a, templates_b, pair_metric; metric_kwargs...)

[MASTER FUNCTION]
Calculates the distance between two multi-molecule configurations using a specified pairwise metric.
It finds the optimal permutation of molecules in `state_a` to match `state_b` by
exhaustively checking all permutations and minimizing the sum of pairwise metric values.
"""
function get_configuration_distance(
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    templates_a::Matrix{Float64},
    templates_b::Matrix{Float64},
    pair_metric::Function;
    metric_kwargs... # To pass things like λ_dist
)
    n_mol = length(state_a)
    @assert n_mol == length(state_b) "States must have the same number of molecules."
    if n_mol < 2
        return 0.0
    end
    
    identity_permutation = 1:n_mol
    pairings = collect(combinations(identity_permutation, 2))
    
    if isempty(pairings)
        return 0.0
    end

    min_sum = Inf
    # Find best permutation of state_a to match state_b by fixing state_b's permutation
    fixed_perm_b = collect(identity_permutation)
    for perm_a in permutations(collect(identity_permutation))
        # Calculate the sum of the metric over all pairs for the current permutation
        current_sum = sum(
            mean( # Note: `mean` assumes the metric returns a vector of per-atom errors.
                pair_metric(
                    templates_a, templates_b,
                    state_a, state_b,
                    [perm_a[p[1]], perm_a[p[2]]], # Get molecules from permuted state_a
                    [fixed_perm_b[p[1]], fixed_perm_b[p[2]]], # Get molecules from fixed state_b
                    ;metric_kwargs...
                )
            ) for p in pairings
        )
        min_sum = min(min_sum, current_sum)
    end

    return min_sum / length(pairings)
end

function get_theta_of_pair(
    template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64},
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    subset_a::Vector{Int}, subset_b::Vector{Int}
)
    # This nested helper correctly uses a specific template with its corresponding state
    function transform_dist(template_centers, state, ss, idx)
        i, j = ss[1], ss[2]
        R1, T1 = state[i]
        R2, T2 = state[j]
        return euclidean(T1 + R1 * template_centers[:,idx], T2 + R2 * template_centers[:,idx])
    end
    
    return [abs(transform_dist(template_centers_a, state_a, subset_a, i) - transform_dist(template_centers_b, state_b, subset_b, i)) for i in 1:size(template_centers_a, 2)]
end

### SCREW AXIS METRIC ###
function get_screw_axis_distance_of_pair(
    templates_a::Matrix{Float64}, templates_b::Matrix{Float64}, # Unused, for signature compatibility
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    subset_a::Vector{Int}, subset_b::Vector{Int};
    λ_dist::Float64 = 0.1
)
    # Use the subsets to extract the specific pair from the full state
    sim_pair_state = [state_a[subset_a[1]], state_a[subset_a[2]]]
    ref_pair_state = [state_b[subset_b[1]], state_b[subset_b[2]]]

    function get_diff_motion(state)
        R1, t1 = state[1]; R2, t2 = state[2]
        return R1' * R2, R1' * (t2 - t1)
    end
    function get_screw_distance(R1, t1, R2, t2)
        diff_R, diff_t = R1' * R2, R1' * (t2 - t1)
        θ = acos(clamp((tr(diff_R) - 1) / 2, -1, 1))
        if isapprox(θ, 0.0, atol=1e-8); return norm(diff_t); end
        axis = [diff_R[3,2] - diff_R[2,3], diff_R[1,3] - diff_R[3,1], diff_R[2,1] - diff_R[1,2]] / (2 * sin(θ))
        d = dot(axis, diff_t)
        return sqrt(λ_dist^2 * d^2 + θ^2)
    end
    
    R_sim, t_sim = get_diff_motion(sim_pair_state)
    R_ref, t_ref = get_diff_motion(ref_pair_state)
    
    # Return as a single-element vector to work with the `mean` call in the master function
    return [get_screw_distance(R_sim, t_sim, R_ref, t_ref)]
end


function get_rmsd_align_one_of_pair(
    templates_a::Matrix{Float64}, templates_b::Matrix{Float64},
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    subset_a::Vector{Int}, subset_b::Vector{Int};
    # No kwargs needed
)
    # Extract the specific pair from the full state using the subsets
    sim_mol_1 = state_a[subset_a[1]]
    sim_mol_2 = state_a[subset_a[2]]
    ref_mol_1 = state_b[subset_b[1]]
    ref_mol_2 = state_b[subset_b[2]]
    
    # --- Helper functions ---
    # function get_coords(state_tuple, templates)
    #     R, T = state_tuple
    #     return [T + R * templates[:, i] for i in 1:size(templates, 2)]
    # end
    function find_superposition_transform(mobile, ref)
        c_mob = mean(mobile); c_ref = mean(ref)
        X = hcat([p .- c_mob for p in mobile]...); Y = hcat([p .- c_ref for p in ref]...)
        H = X * Y'; svd_res = svd(H)
        d = sign(det(svd_res.V * svd_res.U'))
        R = svd_res.V * diagm([1.0, 1.0, d]) * svd_res.U'
        t = c_ref - R * c_mob
        return R, t
    end
    function calculate_superposed_rmsd(mob_A, mob_B, ref_A, ref_B)
        R, t = find_superposition_transform(mob_A, ref_A)
        N = length(mob_A) + length(mob_B)
        sq_err = sum(norm((R * mob_A[i] + t) - ref_A[i])^2 for i in 1:length(ref_A))
        sq_err += sum(norm((R * mob_B[i] + t) - ref_B[i])^2 for i in 1:length(ref_B))
        return sqrt(sq_err / N)
    end
    
    # --- Generate coordinate clouds for the 4 molecules ---
    coords_W = get_point_vector_realization([sim_mol_1], templates_a)
    coords_X = get_point_vector_realization([sim_mol_2], templates_a)
    coords_Y = get_point_vector_realization([ref_mol_1], templates_b)
    coords_Z = get_point_vector_realization([ref_mol_2], templates_b)

    rmsd_A = calculate_superposed_rmsd(coords_W, coords_X, coords_Y, coords_Z)

    # Case 2: Align sim_mol_1 -> ref_mol_2 (and check sim_mol_2 against ref_mol_1)
    rmsd_B = calculate_superposed_rmsd(coords_W, coords_X, coords_Z, coords_Y) # Note swapped ref

    # Return the minimum of the two possibilities for this pair.
    return [min(rmsd_A, rmsd_B)]
end

# These methods simply convert the flat vectors to tuples and call the robust
# core logic from Section 2.

function get_theta_of_pair(
    template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64},
    state_a::Vector{Float64}, state_b::Vector{Float64},
    subset_a::Vector{Int}, subset_b::Vector{Int}
)
    state_a_tuples = convert_flat_state_to_tuples(state_a)
    state_b_tuples = convert_flat_state_to_tuples(state_b)
    return get_theta_of_pair(template_centers_a, template_centers_b, state_a_tuples, state_b_tuples, subset_a, subset_b)
end

function sum_of_permutation(
    template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64},
    state_a::Vector{Float64}, state_b::Vector{Float64},
    perm_a::Vector{Int}, perm_b::Vector{Int},
    all_pairs::Vector{<:AbstractVector}
)
    state_a_tuples = convert_flat_state_to_tuples(state_a)
    state_b_tuples = convert_flat_state_to_tuples(state_b)
    return sum_of_permutation(template_centers_a, template_centers_b, state_a_tuples, state_b_tuples, perm_a, perm_b, all_pairs)
end

function get_configuration_distance(
    template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64},
    state_a::Vector{Float64}, state_b::Vector{Float64}
)
    state_a_tuples = convert_flat_state_to_tuples(state_a)
    state_b_tuples = convert_flat_state_to_tuples(state_b)
    return get_configuration_distance(template_centers_a, template_centers_b, state_a_tuples, state_b_tuples)
end

# Shorthand for same template centers for Vector{Tuple{...}}
function get_configuration_distance(
    template_centers::Matrix{Float64},
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}
)
    return get_configuration_distance(template_centers, template_centers, state_a, state_b)
end

# Shorthand for same template centers for Vector{Float64}
function get_configuration_distance(
    template_centers::Matrix{Float64},
    state_a::Vector{Float64}, state_b::Vector{Float64}
)
    return get_configuration_distance(template_centers, template_centers, state_a, state_b)
end
