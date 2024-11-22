function get_flat_realization(x, template_centers)
    n_mol = length(x) ÷ 6
    [(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]...)...)...]
end

function get_matrix_realization_per_mol(x, template_centers)
    [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]
end

function get_point_realization(x, template_centers)
    n_mol = length(x) ÷ 6
    [Point3f(e) for e in eachcol(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_centers .+ x[i+3:i+5] for i in 1:6:length(x)]...))]
end