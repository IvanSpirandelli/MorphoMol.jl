function get_flat_realization(x, template_mol)
    n_mol = length(x) ÷ 6
    [(hvcat((n_mol), [exp(Rotations.RotationVecGenerator(x[i:i+2]...)) * template_mol .+ x[i+3:i+5] for i in 1:6:length(x)]...)...)...]
end