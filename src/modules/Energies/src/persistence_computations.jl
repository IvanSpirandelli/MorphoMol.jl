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

function get_total_persistence_summed(dim_dgms::Vector{Matrix{Float64}}, weights::Vector{Float64} = [0.1, -0.1, -0.1, 0.0])
    sum([get_persistence(dgm, weight) for (dgm, weight) in zip(dim_dgms, weights)])
end

function get_total_persistence(dgm, weight::Float64 = 1.0)
    weight * sum((dgm[:,2] - dgm[:,1]))
end

function get_divided_persistence(dgm, weight::Float64 = 1.0)
    weight * sum([dgm[i,2] / dgm[i,1] for i in 1:size(dgm)[1]])
end

function get_divided_persistence_summed(dim_dgms::Vector{Matrix{Float64}}, weights::Vector{Float64} = [-0.1, -0.1])
    sum([get_divided_persistence(dgm, weight) for (dgm, weight) in zip(dim_dgms, weights)])
end

function get_interface_diagram_and_filtration(points, n_atoms_per_mol)
    py"""
    import oineus as oin
    import numpy as np
    import torch
    import diode

    def get_interface_diagram_and_filtration(points, n_atoms_per_mol):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        fil = oin.Filtration_double([oin.Simplex_double(s[0], s[1]) for s in simplices])

        def is_multi(sigma):
            return len(set(v // n_atoms_per_mol for v in sigma.vertices)) >= 2
        fil = fil.subfiltration(is_multi)

        def recalculate_filtration_value(cell):
            parts = [v // n_atoms_per_mol for v in cell.vertices]
            n_parts = len(set(parts))
            p_agg = [np.array([0.0, 0.0, 0.0]) for _ in range(n_parts)]
            weights = [0 for _ in range(n_parts)]
            for i, p in enumerate(parts):
                p_agg[p] += points[cell.vertices[i]]
                weights[p] += 1
            bcs = [p_agg[i] / weights[i] for i in range(n_parts)]
            filtration_value = sum([np.linalg.norm(bcs[i] - bcs[j]) for i in range(n_parts) for j in range(i + 1, n_parts)]) / n_parts
            # convert to native Python float from numpy
            return float(filtration_value)

        new_values = [recalculate_filtration_value(cell) for cell in fil]
        fil.set_values(new_values)

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        params.clearing_opt = False
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm, fil
    """
    py"get_interface_diagram_and_filtration"(points, n_atoms_per_mol)
end