using Rotations

function get_flat_realization(x, template_centers)
    if size(template_centers)[2] > 1
        n_mol = length(x) ÷ 6
        [(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]...)...)...]
    else
        return x
    end
end

function get_matrix_realization_per_mol(x, template_centers)
    @assert size(template_centers)[2] > 1
    [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]
end

function get_point3f_realization(x, template_centers)
    if size(template_centers)[2] > 1
        n_mol = length(x) ÷ 6
        return [Point3f(e) for e in eachcol(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]...))]
    else
        n_mol = length(x) ÷ 3
        return [Point3f(e) for e in collect(eachcol(reshape(x, (3,n_mol))))]
    end
end

function get_point_vector_realization(x, template_centers)
    if size(template_centers)[2] > 1
        n_mol = length(x) ÷ 6
        [Vector{Float64}(e) for e in eachcol(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]...))]
    else
        n_mol = length(x) ÷ 3
        return [Vector{Float64}(e) for e in collect(eachcol(reshape(x, (3,n_mol))))]
    end
end