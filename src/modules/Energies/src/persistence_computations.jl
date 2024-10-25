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
    sum([get_total_persistence(dgm, weight) for (dgm, weight) in zip(dim_dgms, weights)])
end

function get_total_persistence(dgm, weight::Float64 = 1.0)
    weight * sum((dgm[:,2] - dgm[:,1]))
end

function get_death_by_birth_persistence(dgm, weight::Float64 = 1.0)
    weight * sum([dgm[i,2] / dgm[i,1] for i in 1:size(dgm)[1]])
end

function get_death_by_birth_persistence_summed(dim_dgms::Vector{Matrix{Float64}}, weights::Vector{Float64} = [-0.1, -0.1])
    sum([get_death_by_birth_persistence(dgm, weight) for (dgm, weight) in zip(dim_dgms, weights)])
end

function get_multichromatic_tetrahedra(points, n_atoms_per_mol)
    py"""
    import numpy as np
    import diode

    def get_multichromatic_tetrahedra(points, n_atoms_per_mol):        
        def is_multi(sigma):
            return len(set(v // n_atoms_per_mol for v in sigma)) >= 2
        points = np.asarray(points)
        tetrahedra = [vs for (vs,fs) in diode.fill_alpha_shapes(points) if len(vs) == 4 and is_multi(vs)]
        return tetrahedra
    """
    py"get_multichromatic_tetrahedra"(points, n_atoms_per_mol)
end

function get_barycenter(points, vertices) 
    Point3f(sum(points[vertices]) / length(vertices))
end

function get_barycentric_subdivision_and_filtration(points, mc_tets)
    barycenters = Vector{Point3f}([])
    filtration = Vector{Tuple{Vector{Int}, Float64}}([])
    total_vertices = 0
    for vs in eachrow(mc_tets)
        part_one = [v+1 for v in vs if div(v, 1206)==0]
        part_two = [v+1 for v in vs if div(v, 1206)==1]
        if length(part_one) == length(part_two) 
            #Get Simplices and Values for even split
            x,y = part_one[1], part_one[2]
            u,v = part_two[1], part_two[2]
            barycenters = [barycenters; [
                get_barycenter(points, [x,y,u,v]),  #1
                get_barycenter(points, [x,v]),      #2
                get_barycenter(points, [x,y,v]),    #3
                get_barycenter(points, [y,v]),      #4
                get_barycenter(points, [y,u,v]),    #5
                get_barycenter(points, [y,u]),      #6
                get_barycenter(points, [x,y,u]),    #7
                get_barycenter(points, [x,u]),      #8
                get_barycenter(points, [x,u,v]),    #9
            ]]
            vertices = [
                ([1] .+ total_vertices, euclidean(get_barycenter(points, [x,y]), get_barycenter(points, [u,v])) / 2.0), 
                ([2] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [v])) / 2.0),
                ([3] .+ total_vertices, euclidean(get_barycenter(points, [x,y]), get_barycenter(points, [v])) / 2.0),
                ([4] .+ total_vertices, euclidean(get_barycenter(points, [y]), get_barycenter(points, [v])) / 2.0),
                ([5] .+ total_vertices, euclidean(get_barycenter(points, [y]), get_barycenter(points, [u, v])) / 2.0),
                ([6] .+ total_vertices, euclidean(get_barycenter(points, [y]), get_barycenter(points, [u])) / 2.0),
                ([7] .+ total_vertices, euclidean(get_barycenter(points, [x,y]), get_barycenter(points, [u])) / 2.0),
                ([8] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [u])) / 2.0),
                ([9] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [v,u])) / 2.0),
            ]
            edges = [
                ([1,2] .+ total_vertices, minimum([vertices[1][2], vertices[2][2]])),
                ([1,3] .+ total_vertices, minimum([vertices[1][2], vertices[3][2]])),
                ([1,4] .+ total_vertices, minimum([vertices[1][2], vertices[4][2]])),
                ([1,5] .+ total_vertices, minimum([vertices[1][2], vertices[5][2]])),
                ([1,6] .+ total_vertices, minimum([vertices[1][2], vertices[6][2]])),
                ([1,7] .+ total_vertices, minimum([vertices[1][2], vertices[7][2]])),
                ([1,8] .+ total_vertices, minimum([vertices[1][2], vertices[8][2]])),
                ([1,9] .+ total_vertices, minimum([vertices[1][2], vertices[9][2]])),
                ([2,3] .+ total_vertices, minimum([vertices[2][2], vertices[3][2]])),
                ([3,4] .+ total_vertices, minimum([vertices[3][2], vertices[4][2]])),
                ([4,5] .+ total_vertices, minimum([vertices[4][2], vertices[5][2]])),
                ([5,6] .+ total_vertices, minimum([vertices[5][2], vertices[6][2]])),
                ([6,7] .+ total_vertices, minimum([vertices[6][2], vertices[7][2]])),
                ([7,8] .+ total_vertices, minimum([vertices[7][2], vertices[8][2]])),
                ([8,9] .+ total_vertices, minimum([vertices[8][2], vertices[9][2]])),
                ([9,2] .+ total_vertices, minimum([vertices[9][2], vertices[2][2]]))
            ]
            triangles = [
                ([1,2,3] .+ total_vertices, minimum([vertices[1][2], vertices[2][2], vertices[3][2]])),
                ([1,3,4] .+ total_vertices, minimum([vertices[1][2], vertices[3][2], vertices[4][2]])),
                ([1,4,5] .+ total_vertices, minimum([vertices[1][2], vertices[4][2], vertices[5][2]])),
                ([1,5,6] .+ total_vertices, minimum([vertices[1][2], vertices[5][2], vertices[6][2]])),
                ([1,6,7] .+ total_vertices, minimum([vertices[1][2], vertices[6][2], vertices[7][2]])),
                ([1,7,8] .+ total_vertices, minimum([vertices[1][2], vertices[7][2], vertices[8][2]])),
                ([1,8,9] .+ total_vertices, minimum([vertices[1][2], vertices[8][2], vertices[9][2]])),
                ([1,9,2] .+ total_vertices, minimum([vertices[1][2], vertices[9][2], vertices[2][2]]))
            ]
            total_vertices += 9
            filtration = [filtration; vertices; edges; triangles]
        end
        if length(part_one) != length(part_two)
            if length(part_one) < length(part_two)
                part_one, part_two = part_two, part_one
            end
            u,v,w = part_one[1], part_one[2], part_one[3]
            x = part_two[1]

            barycenters = [barycenters; [
                get_barycenter(points, [u,v,w,x]),  #1
                get_barycenter(points, [x,v]),      #2
                get_barycenter(points, [x,v,w]),    #3
                get_barycenter(points, [x,w]),      #4
                get_barycenter(points, [w,u,x]),    #5
                get_barycenter(points, [u,x]),      #6
                get_barycenter(points, [u,v,x]),    #7
            ]]
            vertices = [
                ([1] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [u,v,w])) / 2.0), 
                ([2] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [v])) / 2.0),
                ([3] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [v,w])) / 2.0),
                ([4] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [w])) / 2.0),
                ([5] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [w,u])) / 2.0),
                ([6] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [u])) / 2.0),
                ([7] .+ total_vertices, euclidean(get_barycenter(points, [x]), get_barycenter(points, [u,v])) / 2.0),
            ]
            edges = [
                ([1,2] .+ total_vertices, minimum([vertices[1][2], vertices[2][2]])),
                ([1,3] .+ total_vertices, minimum([vertices[1][2], vertices[3][2]])),
                ([1,4] .+ total_vertices, minimum([vertices[1][2], vertices[4][2]])),
                ([1,5] .+ total_vertices, minimum([vertices[1][2], vertices[5][2]])),
                ([1,6] .+ total_vertices, minimum([vertices[1][2], vertices[6][2]])),
                ([1,7] .+ total_vertices, minimum([vertices[1][2], vertices[7][2]])),
                ([2,3] .+ total_vertices, minimum([vertices[2][2], vertices[3][2]])),
                ([3,4] .+ total_vertices, minimum([vertices[3][2], vertices[4][2]])),
                ([4,5] .+ total_vertices, minimum([vertices[4][2], vertices[5][2]])),
                ([5,6] .+ total_vertices, minimum([vertices[5][2], vertices[6][2]])),
                ([6,7] .+ total_vertices, minimum([vertices[6][2], vertices[7][2]])),
                ([7,2] .+ total_vertices, minimum([vertices[7][2], vertices[2][2]]))
            ]
            triangles = [
                ([1,2,3] .+ total_vertices, minimum([vertices[1][2], vertices[2][2], vertices[3][2]])),
                ([1,3,4] .+ total_vertices, minimum([vertices[1][2], vertices[3][2], vertices[4][2]])),
                ([1,4,5] .+ total_vertices, minimum([vertices[1][2], vertices[4][2], vertices[5][2]])),
                ([1,5,6] .+ total_vertices, minimum([vertices[1][2], vertices[5][2], vertices[6][2]])),
                ([1,6,7] .+ total_vertices, minimum([vertices[1][2], vertices[6][2], vertices[7][2]])),
                ([1,7,2] .+ total_vertices, minimum([vertices[1][2], vertices[7][2], vertices[2][2]]))
            ]
            total_vertices += 7
            filtration = [filtration; vertices; edges; triangles]
        end
    end
    barycenters, filtration
end

function calculate_persistence_diagram(upper_star_filtration)
    py"""
    import numpy as np
    import oineus as oin 
    def calculate_persistence_diagram(upper_star_filtration):
        fil = oin.Filtration_double([oin.Simplex_double([v-1 for v in s[0]], s[1]) for s in upper_star_filtration], True)
        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        params.clearing_opt = False
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    py"calculate_persistence_diagram"(upper_star_filtration)
end

function get_interface_with_persistence(points, n_atoms_per_mol)
    mc_tets = get_multichromatic_tetrahedra(points, n_atoms_per_mol)
    barycenters, filtration = get_barycentric_subdivision_and_filtration(points, mc_tets)
    dgms = calculate_persistence_diagram(filtration)
    dgms = [dgms[i] for i in 1:2]
    barycenters, filtration, dgms
end