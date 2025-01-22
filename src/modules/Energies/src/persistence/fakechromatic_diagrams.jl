function get_kic_diagrams(points, subcomplex_labels)
    @assert length(points) >= maximum(subcomplex_labels)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_kic_diagrams(points, subcomplex_labels):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        K = oin.Filtration([oin.Simplex(s, f) for s,f in simplices])

        sls = set(subcomplex_labels)
        sub_complex_simplices = [(s, f) for s,f in simplices if set(s) <= sls]
        L = oin.Filtration([oin.Simplex(s, f) for s,f in sub_complex_simplices])

        kicr = oin.compute_kernel_image_cokernel_reduction(K, L)

        kernel_dgms = kicr.kernel_diagrams()
        image_dgms = kicr.image_diagrams()
        cokernel_dgms = kicr.cokernel_diagrams()
        return kernel_dgms, image_dgms, cokernel_dgms
    """
    # Shifting subcomplex_labels to starting at 0!
    kds, ids, ckds = py"get_kic_diagrams"(points, subcomplex_labels .- 1)
    kernel_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [kds[1], kds[2], kds[3]]]
    image_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [ids[1], ids[2], ids[3]]]
    cokernel_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [ckds[1], ckds[2], ckds[3]]]
    return kernel_dgms, image_dgms, cokernel_dgms
end

function get_kic_diagrams_and_complexes(points, subcomplex_labels)
    @assert length(points) >= maximum(subcomplex_labels)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_kic_diagrams_and_complexes(points, subcomplex_labels):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        K = oin.Filtration([oin.Simplex(s, f) for s,f in simplices])

        sls = set(subcomplex_labels)
        sub_complex_simplices = [(s, f) for s,f in simplices if set(s) <= sls]
        L = oin.Filtration([oin.Simplex(s, f) for s,f in sub_complex_simplices])

        kicr = oin.compute_kernel_image_cokernel_reduction(K, L)

        kernel_dgms = kicr.kernel_diagrams()
        image_dgms = kicr.image_diagrams()
        cokernel_dgms = kicr.cokernel_diagrams()
        return kernel_dgms, image_dgms, cokernel_dgms, K, L
    """    
    # Shifting subcomplex_labels to starting at 0!
    kds, ids, ckds, K, L = py"get_kic_diagrams_and_complexes"(points, subcomplex_labels .- 1)
    kernel_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [kds[1], kds[2], kds[3]]]
    image_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [ids[1], ids[2], ids[3]]]
    cokernel_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [ckds[1], ckds[2], ckds[3]]]
    return kernel_dgms, image_dgms, cokernel_dgms, K, L
end

function get_kic_and_KL_diagrams(points, subcomplex_labels)
    @assert length(points) >= maximum(subcomplex_labels)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_kic_and_KL_diagrams(points, subcomplex_labels):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        K = oin.Filtration([oin.Simplex(s, f) for s,f in simplices])

        sls = set(subcomplex_labels)
        sub_complex_simplices = [(s, f) for s,f in simplices if set(s) <= sls]
        L = oin.Filtration([oin.Simplex(s, f) for s,f in sub_complex_simplices])

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
    # Shifting subcomplex_labels to starting at 0!
    kds, ids, ckds, K_complex_ds, L_complex_ds = py"get_kic_and_KL_diagrams"(points, subcomplex_labels .- 1)
    kernel_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [kds[1], kds[2], kds[3]]]
    image_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [ids[1], ids[2], ids[3]]]
    cokernel_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [ckds[1], ckds[2], ckds[3]]]
    K_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [K_complex_ds[1], K_complex_ds[2], K_complex_ds[3]]]
    L_dgms = [permutedims(hcat([e for e in eachrow(dgm) if !(Inf in e)]...)) for dgm in [L_complex_ds[1], L_complex_ds[2], L_complex_ds[3]]]
    return kernel_dgms, image_dgms, cokernel_dgms, K_dgms, L_dgms
end