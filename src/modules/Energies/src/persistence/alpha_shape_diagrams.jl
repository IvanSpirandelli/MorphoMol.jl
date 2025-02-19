function get_alpha_shape_persistence_diagram(points)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_alpha_shape_persistence_diagram(points):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    asds = py"get_alpha_shape_persistence_diagram"(points)
    [asds[1], asds[2], asds[3]]
end

function debug_alpha_shape(points)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def debug_alpha_shape(points):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])
        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=True)
        return simplices, fil, dgm
    """
    py"debug_alpha_shape"(points)
end