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
        has_a = has_b = False
        for v in sigma.vertices:
            if v < n_atoms_per_mol:
                has_a = True
            else:
                has_b = True
        return has_a and has_b

    fil = fil.subfiltration(is_multi)
    # second argument: True for cohomology, False for homology (incorrect for subfiltration)
    dcmp = oin.Decomposition(fil, True)
    params = oin.ReductionParams()
    params.clearing_opt = False
    dcmp.reduce(params)
    dgm = dcmp.diagram(fil, include_inf_points=False)
    return dgm
"""

function get_interface_diagram(points, n_atoms_per_mol)
    py"get_interface_diagram"(points, n_atoms_per_mol)
end

function get_total_persistence(dim_dgms::Vector{Matrix{Float64}}, weight::Float64 = 1.0)
    sum([get_persistence(dgm, weight) for dgm in dim_dgms])
end

function get_persistence(dgm, weight::Float64 = 1.0)
    weight * sum((dgm[:,2] - dgm[:,1]))
end