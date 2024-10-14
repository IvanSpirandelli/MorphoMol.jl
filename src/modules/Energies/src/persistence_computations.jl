function get_interface_diagram(points, n_atoms_per_mol)
    py"""
    import oineus as oin
    import numpy as np
    import torch
    import diode

    def get_interface_diagram(points, n_atoms_per_mol):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        fil = oin.Filtration_double([oin.Simplex_double(s[0], s[1]) for s in simplices])

        def is_multi(sigma):
            return len(set(v // n_atoms_per_mol for v in sigma.vertices)) >= 2

        fil = fil.subfiltration(is_multi)
        # second argument: True for cohomology, False for homology (incorrect for subfiltration)
        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        params.clearing_opt = False
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    py"get_interface_diagram"(points, n_atoms_per_mol)
end

function get_persistence_diagram(points)
    py"""
    import oineus as oin
    import numpy as np
    import torch
    import diode

    def get_persistence_diagram(points):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        fil = oin.Filtration_double([oin.Simplex_double(s[0], s[1]) for s in simplices])

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    py"get_persistence_diagram"(points)
end

function get_total_persistence(dim_dgms::Vector{Matrix{Float64}}, weights::Vector{Float64} = [0.1, -0.1, -0.1, 0.0])
    sum([get_persistence(dgm, weight) for (dgm, weight) in zip(dim_dgms, weights)])
end

function get_persistence(dgm, weight::Float64 = 1.0)
    weight * sum((dgm[:,2] - dgm[:,1]))
end