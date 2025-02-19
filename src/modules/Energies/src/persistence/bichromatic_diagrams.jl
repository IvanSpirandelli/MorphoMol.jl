function bichromatic_persistence_all_diagrams(points, colors, sub_complex_color)
    py"""
    import math
    import numpy as np
    import oineus as oin
    import bichromatic_delaunay as bd
    def bichromatic_persistence_all_diagrams(points, colors, sub_complex_color):
        simplices, values = bd.bichromatic_delaunay(points, colors)
        K = oin.Filtration([oin.Simplex(s, fv) for s,fv in zip(simplices, values)])
        L = oin.Filtration([oin.Simplex(s, fv) for s,fv in zip(simplices, values) if all(math.floor(colors[e]) == sub_complex_color for e in s)])

        kicr = oin.compute_kernel_image_cokernel_reduction(K, L)

        kernel_dgms = kicr.kernel_diagrams()
        image_dgms = kicr.image_diagrams()
        cokernel_dgms = kicr.cokernel_diagrams()

        params = oin.ReductionParams()
        params.clearing_opt = False

        dcmp_K = oin.Decomposition(K, False) #True means dualize, i.e. cohomology
        dcmp_K.reduce(params)
        K_dgms = dcmp_K.diagram(K, include_inf_points=False)

        dcmp_L = oin.Decomposition(L, False) 
        dcmp_L.reduce(params)
        L_dgms = dcmp_L.diagram(L, include_inf_points=False)

        return kernel_dgms, image_dgms, cokernel_dgms, K_dgms, L_dgms
    """
    kds, ids, ckds, K_complex_ds, L_complex_ds = py"bichromatic_persistence_all_diagrams"(points, colors, sub_complex_color)
    kernel_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [kds[1], kds[2], kds[3]]]
    image_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [ids[1], ids[2], ids[3]]]
    cokernel_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [ckds[1], ckds[2], ckds[3]]]
    K_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [K_complex_ds[1], K_complex_ds[2], K_complex_ds[3]]]
    L_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [L_complex_ds[1], L_complex_ds[2], L_complex_ds[3]]]
    return kernel_dgms, image_dgms, cokernel_dgms, K_dgms, L_dgms
end

function bichromatic_delaunay(points, colors)
    py"""
    import numpy as np
    import oineus as oin
    import bichromatic_delaunay as bd
    def bichromatic_persistence_all_diagrams(points, colors):
        simplices, values = bd.bichromatic_delaunay(points, colors)
        return simplices, values
    """
    py"bichromatic_persistence_all_diagrams"(points, colors)
end
    