function get_matched_distances_between_transformation_offsets(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64})
    function transform_dist(template_centers, state, idx)
        i = 1
        j = 2
        R1 = exp(RotationVecGenerator(state[(i-1) * 6 + 1:(i-1) * 6 + 3]...))
        T1 = state[(i-1) * 6 + 4:(i-1) * 6 + 6]
        R2 = exp(RotationVecGenerator(state[(j-1) * 6 + 1:(j-1) * 6 + 3]...))
        T2 = state[(j-1) * 6 + 4:(j-1) * 6 + 6]
        euclidean(T1 + R1 * template_centers[:,idx], T2 + R2 * template_centers[:,idx])
    end

    [abs(transform_dist(template_centers_a, state_a, i) - transform_dist(template_centers_b, state_b, i)) for i in 1:size(template_centers_a)[2]]
end

function average_offset_distance(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64})
    n = size(template_centers_a)[2]
    sum(get_matched_distances_between_transformation_offsets(template_centers_a, template_centers_b, state_a, state_b)) / n
end

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

function sum_of_permutation(template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64}, state_a::Vector{Float64}, state_b::Vector{Float64}, perm_a::Vector{Int}, perm_b::Vector{Int})
    n = size(template_centers_a)[2]
    d1 = sum(get_matched_distances_between_transformation_offsets(template_centers_a, template_centers_b, state_a, state_b, perm_a[[1,2]], perm_b[[1,2]])) / n
    d2 = sum(get_matched_distances_between_transformation_offsets(template_centers_a, template_centers_b, state_a, state_b, perm_a[[2,3]], perm_b[[2,3]])) / n
    d3 = sum(get_matched_distances_between_transformation_offsets(template_centers_a, template_centers_b, state_a, state_b, perm_a[[1,3]], perm_b[[1,3]])) / n
    println("$(d1), $(d2), $(d3)")
    (d1 + d2 + d3)/3.0
end