function solvation_free_energy(
    prefactors::AbstractVector, 
    atom_coordinates::Vector, 
    molecule_size::Int, 
    atom_radii::Vector, 
    probe_radius::Float64, 
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64,
    delaunay_eps::Float64 = 1.0)
    measures = get_hadwiger_measures_and_linear_overlap(
        atom_coordinates,
        molecule_size,
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
    sum(measures .* prefactors)
end
