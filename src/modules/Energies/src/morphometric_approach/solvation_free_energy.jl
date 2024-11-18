function solvation_free_energy(
    atom_coordinates::Vector, 
    atom_radii::Vector, 
    probe_radius::Float64, 
    prefactors::AbstractVector, 
    delaunay_eps::Float64 = 1.0)
    measures = get_geometric_measures(
        atom_coordinates,
        atom_radii,
        probe_radius,
        delaunay_eps
    )
    sum(measures .* prefactors)
end

function solvation_free_energy(
    atom_coordinates::Vector, 
    molecule_size::Int, 
    atom_radii::Vector, 
    probe_radius::Float64, 
    prefactors::AbstractVector, 
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64,
    delaunay_eps::Float64 = 1.0)
    measures = get_geometric_measures_and_overlap_value(
        atom_coordinates,
        molecule_size,
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
    sum(measures .* [prefactors; 1.0])
end