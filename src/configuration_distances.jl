using Combinatorics

# General case:
function get_matched_distances_between_transformation_offsets(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, subset_a::Vector{Int}, subset_b::Vector{Int})
    function transform_dist(template_centers, state, ss, idx)
        i = ss[1]
        j = ss[2]
        R1 = exp(RotationVecGenerator(state[(i-1) * 6 + 1:(i-1) * 6 + 3]...))
        T1 = state[(i-1) * 6 + 4:(i-1) * 6 + 6]
        R2 = exp(RotationVecGenerator(state[(j-1) * 6 + 1:(j-1) * 6 + 3]...))
        T2 = state[(j-1) * 6 + 4:(j-1) * 6 + 6]
        euclidean(T1 + R1 * template_centers[:,idx], T2 + R2 * template_centers[:,idx])
    end
    [abs(transform_dist(template_centers_a, state_a, subset_a, i) - transform_dist(template_centers_b, state_b, subset_b, i)) for i in 1:size(template_centers_a)[2]]
end

function sum_of_permutation(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, perm_a::Vector{Int}, perm_b::Vector{Int}, all_pairs::Vector{Vector{Int}}, n_mol::Int)
    n = size(template_centers_a)[2]
    matched_distances_by_pair = [get_matched_distances_between_transformation_offsets(template_centers_a, template_centers_b, state_a, state_b, perm_a[p], perm_b[p]) for p in all_pairs]
    spreads_by_by_pair = [maximum(mdbp) - minimum(mdbp) for mdbp in matched_distances_by_pair]
    #println([(maximum(mdbp), minimum(mdbp)) for mdbp in matched_distances_by_pair])
    sum(sum(mdbp)/ n for mdbp in matched_distances_by_pair) / binomial(n_mol,2), sum(spreads_by_by_pair) / binomial(n_mol,2)
end

function average_offset_distance(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, n_mol::Int)
    identity_permutation = collect(1:n_mol)
    pairings = collect(combinations(identity_permutation, 2))
    sum_of_matched_distances_and_spreads = [sum_of_permutation(template_centers_a, template_centers_b, state_a, state_b, perm, identity_permutation, pairings, n_mol) for perm in collect(permutations(identity_permutation))]
    min_index = argmin([x[1] for x in sum_of_matched_distances_and_spreads])
    sum_of_matched_distances_and_spreads[min_index]
end

function get_matched_distances_between_transformation_offsets(template_centers::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, subset_a::Vector{Int}, subset_b::Vector{Int})
    get_matched_distances_between_transformation_offsets(template_centers, template_centers, state_a, state_b, subset_a, subset_b)
end

function sum_of_permutation(template_centers::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, perm_a::Vector{Int}, perm_b::Vector{Int}, all_pairs::Vector{Vector{Int}}, n_mol::Int)
    sum_of_permutation(template_centers, template_centers, state_a, state_b, perm_a, perm_b, all_pairs, n_mol)
end

function average_offset_distance(template_centers::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, n_mol::Int)
    average_offset_distance(template_centers, template_centers, state_a, state_b, n_mol)
end

# function get_rotation_distance(R1, R2)
#     diffR = R1 * R2'
#     acos((tr(diffR)-1.0) * 0.5)
# end

# function get_matched_rotation_offset_distance(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, subset_a::Vector{Int}, subset_b::Vector{Int})
#     function get_diff_rotation(template_centers, state, ss)
#         i = ss[1]
#         j = ss[2]
#         R1 = exp(RotationVecGenerator(state[(i-1) * 6 + 1:(i-1) * 6 + 3]...))
#         R2 = exp(RotationVecGenerator(state[(j-1) * 6 + 1:(j-1) * 6 + 3]...))
#         R1 * R2'
#     end
#     abs(get_rotation_distance(get_diff_rotation(template_centers_a, state_a, subset_a), get_diff_rotation(template_centers_b, state_b, subset_b)))
# end

# function sum_of_rotation_offsets_of_permutation(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, perm_a::Vector{Int}, perm_b::Vector{Int}, all_pairs::Vector{Vector{Int}}, n_mol::Int)
#     sum(get_matched_rotation_offset_distance(template_centers_a, template_centers_b, state_a, state_b, perm_a[p], perm_b[p]) for p in all_pairs) / binomial(n_mol,2)
# end

# function average_rotation_offset_distance(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, n_mol::Int)
#     @assert n_mol == 2 "Only tested for n_mol == 2"
#     identity_permutation = collect(1:n_mol)
#     pairings = collect(combinations(identity_permutation, 2))
#     minimum([sum_of_rotation_offsets_of_permutation(template_centers_a, template_centers_b, state_a, state_b, perm, identity_permutation, pairings, n_mol) for perm in collect(permutations(identity_permutation))])
# end