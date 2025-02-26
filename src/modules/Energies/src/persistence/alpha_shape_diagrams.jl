function get_alpha_shape_persistence_diagram(points, exact = false)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_alpha_shape_persistence_diagram(points, exact):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points, exact=exact)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    asds = py"get_alpha_shape_persistence_diagram"(points, exact)
    [asds[1], asds[2], asds[3]]
end

function debug_alpha_shape(points)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def debug_alpha_shape(points):
        points = np.asarray(points)
        np.save('case2.npy', points)
        simplices = diode.fill_alpha_shapes(points)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])
        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return simplices, fil, dgm
    """
    py"debug_alpha_shape"(points)
end