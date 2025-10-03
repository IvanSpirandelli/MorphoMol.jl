using Rotations

function get_flat_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    n_mol = length(x)
    [(hvcat((n_mol), [R * tc .+ t for ((R,t), tc) in zip(x, template_centers)]...)...)...]
end

function get_point_vector_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    n_mol = length(x)
    [Vector{Float64}(e) for e in eachcol(hvcat((n_mol), [R * tc .+ t for ((R,t), tc) in zip(x, template_centers)]...))]
end

function get_flat_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64})
    if size(template_centers)[2] > 1
        n_mol = length(x)
        [(hvcat((n_mol), [R * template_centers .+ t for (R,t) in x]...)...)...]
    else
        @assert false "This function is not implemented for single hard spheres."
    end
end

function get_matrix_realization_per_mol(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64})
    @assert size(template_centers)[2] > 1
    [R * template_centers .+ t for (R,t) in x]
end

function get_point_vector_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64})
    if size(template_centers)[2] > 1
        n_mol = length(x)
        [Vector{Float64}(e) for e in eachcol(hvcat((n_mol), [R * template_centers .+ t for (R,t) in x]...))]
    else
        @assert false "This function is not implemented for single hard spheres."
    end
end

function get_point3f_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64})
    if size(template_centers)[2] > 1
        n_mol = length(x)
        return [Point3f(e) for e in eachcol(hvcat((n_mol), [R * template_centers .+ t for (R,t) in x]...))]
    else
        @assert false "This function is not implemented for single hard spheres."
    end
end


# These implemenations sampled Rotations badly. Kept them for evaluation of old data.
function get_flat_realization(x::Vector{Float64}, template_centers::Matrix{Float64})
    if size(template_centers)[2] > 1
        n_mol = length(x) ÷ 6
        [(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]...)...)...]
    else
        return x
    end
end

function get_matrix_realization_per_mol(x::Vector{Float64}, template_centers::Matrix{Float64})
    @assert size(template_centers)[2] > 1
    [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]
end

function get_point3f_realization(x::Vector{Float64}, template_centers::Matrix{Float64})
    if size(template_centers)[2] > 1
        n_mol = length(x) ÷ 6
        return [Point3f(e) for e in eachcol(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]...))]
    else
        n_mol = length(x) ÷ 3
        return [Point3f(e) for e in collect(eachcol(reshape(x, (3,n_mol))))]
    end
end

function get_point_vector_realization(x::Vector{Float64}, template_centers::Matrix{Float64})
    if size(template_centers)[2] > 1
        n_mol = length(x) ÷ 6
        [Vector{Float64}(e) for e in eachcol(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]...))]
    else
        n_mol = length(x) ÷ 3
        return [Vector{Float64}(e) for e in collect(eachcol(reshape(x, (3,n_mol))))]
    end
end