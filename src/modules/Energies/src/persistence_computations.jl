function get_alpha_shape_persistence_diagram(points)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_alpha_shape_persistence_diagram(points):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points)
        fil = oin.Filtration_double([oin.Simplex_double(s[0], s[1]) for s in simplices])

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    py"get_alpha_shape_persistence_diagram"(points)
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

function get_interface_persistence_diagram_from_upper_star_filtration(filtration)
    py"""
    import numpy as np
    import oineus as oin 
    def get_interface_persistence_diagram_from_upper_star_filtration(filtration):
        fil = oin.Filtration_double([oin.Simplex_double([v-1 for v in s[0]], s[1]) for s in filtration], True) #True means negate, i.e. upper/lower star???
        dcmp = oin.Decomposition(fil, False) #True means dualize, i.e. cohomology
        params = oin.ReductionParams()
        params.clearing_opt = False
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    py"get_interface_persistence_diagram_from_upper_star_filtration"(filtration)
end

function get_interface_persistence_diagram_and_geometry(points, n_atoms_per_mol)
    mc_tets = get_multichromatic_tetrahedra(points, n_atoms_per_mol)
    barycenters, filtration = get_barycentric_subdivision_and_filtration(points, mc_tets, n_atoms_per_mol)
    dgms = get_interface_persistence_diagram_from_upper_star_filtration(filtration)
    dgms = [dgms[i] for i in 1:2]
    dgms, barycenters, filtration
end

function get_interface_persistence_diagram(points, n_atoms_per_mol)
    mc_tets = get_multichromatic_tetrahedra(points, n_atoms_per_mol)
    _, filtration = get_barycentric_subdivision_and_filtration(points, mc_tets, n_atoms_per_mol)
    dgms = get_interface_persistence_diagram_from_upper_star_filtration(filtration)
    dgms = [dgms[i] for i in 1:2]
    dgms
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
    tets = py"get_multichromatic_tetrahedra"(points, n_atoms_per_mol)
    tets = [t .+ 1 for t in tets]
end

function get_chromatic_partitioning(tet, n_atoms_per_mol)
    divs = [div(v - 1, n_atoms_per_mol) for v in tet]
    parts = Dict(k => Vector{Int}([]) for k in unique(divs))
    for (i, d) in enumerate(divs)
        push!(parts[d], tet[i])
    end
    sort([v for v in values(parts)], by=length, rev = true)
end

function get_barycenter(points, vertices) 
    Point3f(sum(points[vertices]) / length(vertices))
end

function get_or_create_bc_simplex_id_and_val!(partitioning::Vector{Vector{Int}}, uob_to_barycenter_simplices, points)
    simplex = Tuple(sort(vcat(partitioning...)))
    try
        (id, val) = uob_to_barycenter_simplices[simplex]
        return true, (id, val)
    catch
        id = length(uob_to_barycenter_simplices) + 1
        if length(partitioning) == 2
            part_one, part_two = partitioning
            val = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_two)) / 2.0
        elseif length(partitioning) == 3
            part_one, part_two, part_three = partitioning
            a = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_two))
            b = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_three))
            c = euclidean(get_barycenter(points, part_two), get_barycenter(points, part_three))
            val = (a + b + c) / 3.0
        elseif length(partitioning) == 4
            part_one, part_two, part_three, part_four = partitioning
            a = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_two))
            b = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_three))
            c = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_four))
            d = euclidean(get_barycenter(points, part_two), get_barycenter(points, part_three))
            e = euclidean(get_barycenter(points, part_two), get_barycenter(points, part_four))
            f = euclidean(get_barycenter(points, part_three), get_barycenter(points, part_four))
            val = (a + b + c + d + e + f) / 6.0
        end
        uob_to_barycenter_simplices[simplex] = (id, val)
        return false, (id, val)
    end
end 

function get_barycentric_subdivision_and_filtration(points, mc_tets, n_atoms_per_mol::Int)
    barycenters = Vector{Point3f}([])
    filtration = Set{Tuple{Vector{Int32}, Float64}}([])
    uob_to_barycenter_simplices = Dict{Any, Any}()
    for vs in eachrow(mc_tets)
        parts = get_chromatic_partitioning(vs, n_atoms_per_mol)
        if length(parts) == 2
            part_one, part_two = parts
            if length(part_one) == length(part_two) == 2
                #Get Simplices and Values for even split
                u, v = part_one
                x, y = part_two
    
                created = Vector{Bool}([])
                vertices = Vector{Tuple{Int, Float64}}([])
                for uob_simplex in [[[u,v],[x,y]], [[x],[v]], [[x,y],[v]], [[y],[v]], [[y],[u,v]], [[y],[u]], [[x,y],[u]], [[x],[u]], [[x],[u,v]]]
                    exists, (id, val) = get_or_create_bc_simplex_id_and_val!(uob_simplex, uob_to_barycenter_simplices, points)
                    push!(vertices, (id, val))
                    push!(created, !exists)
                end
    
                #Placing the points for visualization
                bc_e = [
                    get_barycenter(points, [x,v]),
                    get_barycenter(points, [y,v]),
                    get_barycenter(points, [y,u]),
                    get_barycenter(points, [x,u])
                ]
    
                barycenters = [barycenters; [
                    get_barycenter(bc_e, [1,2,3,4]),    #u v x y
                    bc_e[1],                            #x v
                    get_barycenter(bc_e, [1,2]),        #x y v
                    bc_e[2],                            #y v
                    get_barycenter(bc_e, [2,3]),        #y u v
                    bc_e[3],                            #y u
                    get_barycenter(bc_e, [3,4]),        #x y u
                    bc_e[4],                            #x u
                    get_barycenter(bc_e, [1,4]),        #x u v
                ][created]]
    
                # We sort edges to filter out duplicates, there will never be dubplicate triangles and vertices dont need to be sorted.
                edges = [
                    (sort!([vertices[1][1], vertices[2][1]]), minimum([vertices[1][2], vertices[2][2]])),
                    (sort!([vertices[1][1], vertices[3][1]]), minimum([vertices[1][2], vertices[3][2]])),
                    (sort!([vertices[1][1], vertices[4][1]]), minimum([vertices[1][2], vertices[4][2]])),
                    (sort!([vertices[1][1], vertices[5][1]]), minimum([vertices[1][2], vertices[5][2]])),
                    (sort!([vertices[1][1], vertices[6][1]]), minimum([vertices[1][2], vertices[6][2]])),
                    (sort!([vertices[1][1], vertices[7][1]]), minimum([vertices[1][2], vertices[7][2]])),
                    (sort!([vertices[1][1], vertices[8][1]]), minimum([vertices[1][2], vertices[8][2]])),
                    (sort!([vertices[1][1], vertices[9][1]]), minimum([vertices[1][2], vertices[9][2]])),
                    (sort!([vertices[2][1], vertices[3][1]]), minimum([vertices[2][2], vertices[3][2]])),
                    (sort!([vertices[3][1], vertices[4][1]]), minimum([vertices[3][2], vertices[4][2]])),
                    (sort!([vertices[4][1], vertices[5][1]]), minimum([vertices[4][2], vertices[5][2]])),
                    (sort!([vertices[5][1], vertices[6][1]]), minimum([vertices[5][2], vertices[6][2]])),
                    (sort!([vertices[6][1], vertices[7][1]]), minimum([vertices[6][2], vertices[7][2]])),
                    (sort!([vertices[7][1], vertices[8][1]]), minimum([vertices[7][2], vertices[8][2]])),
                    (sort!([vertices[8][1], vertices[9][1]]), minimum([vertices[8][2], vertices[9][2]])),
                    (sort!([vertices[9][1], vertices[2][1]]), minimum([vertices[9][2], vertices[2][2]]))
                ]
    
                triangles = [
                    ([vertices[1][1], vertices[2][1], vertices[3][1]], minimum([vertices[1][2], vertices[2][2], vertices[3][2]])),
                    ([vertices[1][1], vertices[3][1], vertices[4][1]], minimum([vertices[1][2], vertices[3][2], vertices[4][2]])),
                    ([vertices[1][1], vertices[4][1], vertices[5][1]], minimum([vertices[1][2], vertices[4][2], vertices[5][2]])),
                    ([vertices[1][1], vertices[5][1], vertices[6][1]], minimum([vertices[1][2], vertices[5][2], vertices[6][2]])),
                    ([vertices[1][1], vertices[6][1], vertices[7][1]], minimum([vertices[1][2], vertices[6][2], vertices[7][2]])),
                    ([vertices[1][1], vertices[7][1], vertices[8][1]], minimum([vertices[1][2], vertices[7][2], vertices[8][2]])),
                    ([vertices[1][1], vertices[8][1], vertices[9][1]], minimum([vertices[1][2], vertices[8][2], vertices[9][2]])),
                    ([vertices[1][1], vertices[9][1], vertices[2][1]], minimum([vertices[1][2], vertices[9][2], vertices[2][2]]))
                ]
                union!(filtration, Set([([v[1]], v[2]) for v in vertices]))
                union!(filtration, Set(edges))
                union!(filtration, Set(triangles))
            elseif 3 == length(part_one) != length(part_two) == 1
                u,v,w = part_one[1], part_one[2], part_one[3]
                x = part_two[1]
    
                created = Vector{Bool}([])
                vertices = Vector{Tuple{Int, Float64}}([])
                for uob_simplex in [[[u,v,w],[x]], [[x],[v]], [[x],[w,v]], [[w],[x]], [[w,u],[x]], [[u],[x]], [[u,v],[x]]]
                    exists, (id, val) = get_or_create_bc_simplex_id_and_val!(uob_simplex, uob_to_barycenter_simplices, points)
                    push!(vertices, (id, val))
                    push!(created, !exists)
                end
    
                bc_e = [
                    get_barycenter(points, [x,v]),      
                    get_barycenter(points, [x,w]),      
                    get_barycenter(points, [u,x])       
                ]
    
                barycenters = [barycenters; [
                    get_barycenter(bc_e, [1,2,3]),  #u v w x
                    bc_e[1],                        #v x
                    get_barycenter(bc_e, [1,2]),    #v w x
                    bc_e[2],                        #w x
                    get_barycenter(bc_e, [2,3]),    #w u x
                    bc_e[3],                        #u x
                    get_barycenter(bc_e, [1,3]),    #u v x
                ][created]]
                
                # We sort edges to filter out duplicates, there will never be dubplicate triangles and vertices dont need to be sorted.
                edges = [
                    (sort!([vertices[1][1], vertices[2][1]]), minimum([vertices[1][2], vertices[2][2]])),
                    (sort!([vertices[1][1], vertices[3][1]]), minimum([vertices[1][2], vertices[3][2]])),
                    (sort!([vertices[1][1], vertices[4][1]]), minimum([vertices[1][2], vertices[4][2]])),
                    (sort!([vertices[1][1], vertices[5][1]]), minimum([vertices[1][2], vertices[5][2]])),
                    (sort!([vertices[1][1], vertices[6][1]]), minimum([vertices[1][2], vertices[6][2]])),
                    (sort!([vertices[1][1], vertices[7][1]]), minimum([vertices[1][2], vertices[7][2]])),
                    (sort!([vertices[2][1], vertices[3][1]]), minimum([vertices[2][2], vertices[3][2]])),
                    (sort!([vertices[3][1], vertices[4][1]]), minimum([vertices[3][2], vertices[4][2]])),
                    (sort!([vertices[4][1], vertices[5][1]]), minimum([vertices[4][2], vertices[5][2]])),
                    (sort!([vertices[5][1], vertices[6][1]]), minimum([vertices[5][2], vertices[6][2]])),
                    (sort!([vertices[6][1], vertices[7][1]]), minimum([vertices[6][2], vertices[7][2]])),
                    (sort!([vertices[7][1], vertices[2][1]]), minimum([vertices[7][2], vertices[2][2]]))
                ]
    
                triangles = [
                    ([vertices[1][1], vertices[2][1], vertices[3][1]], minimum([vertices[1][2], vertices[2][2], vertices[3][2]])),
                    ([vertices[1][1], vertices[3][1], vertices[4][1]], minimum([vertices[1][2], vertices[3][2], vertices[4][2]])),
                    ([vertices[1][1], vertices[4][1], vertices[5][1]], minimum([vertices[1][2], vertices[4][2], vertices[5][2]])),
                    ([vertices[1][1], vertices[5][1], vertices[6][1]], minimum([vertices[1][2], vertices[5][2], vertices[6][2]])),
                    ([vertices[1][1], vertices[6][1], vertices[7][1]], minimum([vertices[1][2], vertices[6][2], vertices[7][2]])),
                    ([vertices[1][1], vertices[7][1], vertices[2][1]], minimum([vertices[1][2], vertices[7][2], vertices[2][2]]))
                ]
    
                union!(filtration, Set([([v[1]], v[2]) for v in vertices]))
                union!(filtration, Set(edges))
                union!(filtration, Set(triangles))
            else
                error("Invalid partitioning: $vs -> $parts")
            end
        elseif length(parts) == 3
            part_one, part_two, part_three = parts
            a, b = part_one
            u = part_two[1]
            x = part_three[1]

            created = Vector{Bool}([])
            vertices = Vector{Tuple{Int, Float64}}([])
            for uob_simplex in [[[a,b],[u],[x]], [[a],[u]], [[a],[u],[x]], [[a],[x]], [[a,b],[x]], [[b],[x]], [[b],[u],[x]], [[b],[u]], [[a,b],[u]], [[u],[x]]]
                exists, (id, val) = get_or_create_bc_simplex_id_and_val!(uob_simplex, uob_to_barycenter_simplices, points)
                push!(vertices, (id, val))
                push!(created, !exists)
            end

            bc_e = [
                get_barycenter(points, [a,u]), #1
                get_barycenter(points, [a,x]), #2   
                get_barycenter(points, [b,x]), #3
                get_barycenter(points, [b,u]), #4
                get_barycenter(points, [u,x])  #5
            ]

            barycenters = [barycenters; [
                get_barycenter(bc_e, [1,2,3,4,5]),  #a b u x
                bc_e[1],                            #a u 
                get_barycenter(bc_e, [1,2,5]),      #a u x
                bc_e[2],                            #a x
                get_barycenter(bc_e, [2,3]),        #a b x
                bc_e[3],                            #b x
                get_barycenter(bc_e, [3,4,5]),      #b u x
                bc_e[4],                            #b u
                get_barycenter(bc_e, [1,4]),        #a b u
                bc_e[5],                            #u x
            ][created]]

            # We sort edges to filter out duplicates, there will never be dubplicate triangles and vertices dont need to be sorted.
            edges = [
                (sort!([vertices[1][1], vertices[2][1]]), minimum([vertices[1][2], vertices[2][2]])),
                (sort!([vertices[1][1], vertices[3][1]]), minimum([vertices[1][2], vertices[3][2]])),
                (sort!([vertices[1][1], vertices[4][1]]), minimum([vertices[1][2], vertices[4][2]])),
                (sort!([vertices[1][1], vertices[5][1]]), minimum([vertices[1][2], vertices[5][2]])),
                (sort!([vertices[1][1], vertices[6][1]]), minimum([vertices[1][2], vertices[6][2]])),
                (sort!([vertices[1][1], vertices[7][1]]), minimum([vertices[1][2], vertices[7][2]])),
                (sort!([vertices[1][1], vertices[8][1]]), minimum([vertices[1][2], vertices[8][2]])),
                (sort!([vertices[1][1], vertices[9][1]]), minimum([vertices[1][2], vertices[9][2]])),
                (sort!([vertices[2][1], vertices[3][1]]), minimum([vertices[2][2], vertices[3][2]])),
                (sort!([vertices[3][1], vertices[4][1]]), minimum([vertices[3][2], vertices[4][2]])),
                (sort!([vertices[4][1], vertices[5][1]]), minimum([vertices[4][2], vertices[5][2]])),
                (sort!([vertices[5][1], vertices[6][1]]), minimum([vertices[5][2], vertices[6][2]])),
                (sort!([vertices[6][1], vertices[7][1]]), minimum([vertices[6][2], vertices[7][2]])),
                (sort!([vertices[7][1], vertices[8][1]]), minimum([vertices[7][2], vertices[8][2]])),
                (sort!([vertices[8][1], vertices[9][1]]), minimum([vertices[8][2], vertices[9][2]])),
                (sort!([vertices[9][1], vertices[2][1]]), minimum([vertices[9][2], vertices[2][2]])),
                (sort!([vertices[1][1], vertices[10][1]]), minimum([vertices[1][2], vertices[10][2]])),
                (sort!([vertices[3][1], vertices[10][1]]), minimum([vertices[3][2], vertices[10][2]])),
                (sort!([vertices[7][1], vertices[10][1]]), minimum([vertices[7][2], vertices[10][2]])),
            ]

            triangles = [
                ([vertices[1][1], vertices[2][1], vertices[3][1]], minimum([vertices[1][2], vertices[2][2], vertices[3][2]])),
                ([vertices[1][1], vertices[3][1], vertices[4][1]], minimum([vertices[1][2], vertices[3][2], vertices[4][2]])),
                ([vertices[1][1], vertices[4][1], vertices[5][1]], minimum([vertices[1][2], vertices[4][2], vertices[5][2]])),
                ([vertices[1][1], vertices[5][1], vertices[6][1]], minimum([vertices[1][2], vertices[5][2], vertices[6][2]])),
                ([vertices[1][1], vertices[6][1], vertices[7][1]], minimum([vertices[1][2], vertices[6][2], vertices[7][2]])),
                ([vertices[1][1], vertices[7][1], vertices[8][1]], minimum([vertices[1][2], vertices[7][2], vertices[8][2]])),
                ([vertices[1][1], vertices[8][1], vertices[9][1]], minimum([vertices[1][2], vertices[8][2], vertices[9][2]])),
                ([vertices[1][1], vertices[9][1], vertices[2][1]], minimum([vertices[1][2], vertices[9][2], vertices[2][2]])),
                ([vertices[1][1], vertices[3][1], vertices[10][1]], minimum([vertices[1][2], vertices[3][2], vertices[10][2]])),
                ([vertices[1][1], vertices[7][1], vertices[10][1]], minimum([vertices[1][2], vertices[7][2], vertices[10][2]])),
            ]
            union!(filtration, Set([([v[1]], v[2]) for v in vertices]))
            union!(filtration, Set(edges))
            union!(filtration, Set(triangles))
        elseif length(parts) == 4
            part_one, part_two, part_three, part_four = parts
            a = part_one[1]
            i = part_two[1]
            u = part_three[1]
            x = part_four[1]

            created = Vector{Bool}([])
            vertices = Vector{Tuple{Int, Float64}}([])
            for uob_simplex in [[[a],[i]], [[a],[u]], [[a],[x]], [[i],[u]], [[i],[x]], [[u],[x]], [[a],[i],[u]], [[a],[i],[x]], [[i],[u],[x]],  [[a],[u],[x]], [[a],[i],[u],[x]]]
                exists, (id, val) = get_or_create_bc_simplex_id_and_val!(uob_simplex, uob_to_barycenter_simplices, points)
                push!(vertices, (id, val))
                push!(created, !exists)
            end

            bc_e = [
                get_barycenter(points, [a,i]), #1
                get_barycenter(points, [a,u]), #2   
                get_barycenter(points, [a,x]), #3
                get_barycenter(points, [i,u]), #4
                get_barycenter(points, [i,x]), #5
                get_barycenter(points, [u,x]), #6
            ]

            barycenters = [barycenters; [
                bc_e[1],                                #a i            1
                bc_e[2],                                #a u            2
                bc_e[3],                                #a x            3
                bc_e[4],                                #i u            4
                bc_e[5],                                #i x            5
                bc_e[6],                                #u x            6
                get_barycenter(bc_e, [1,2,4]),          #a i u          7
                get_barycenter(bc_e, [1,3,5]),          #a i x          8
                get_barycenter(bc_e, [4,5,6]),          #i u x          9
                get_barycenter(bc_e, [2,3,6]),          #a u x          10  
                get_barycenter(bc_e, [1,2,3,4,5,6]),    #a i u x        11
            ][created]]

            # We sort edges to filter out duplicates, there will never be dubplicate triangles and vertices dont need to be sorted.
            edges = Vector{Tuple{Vector{Int}, Float64}}([])
            for (i, j) in [(1,7), (1,8), (1,11), (2,7), (2,10), (2,11), (3,8), (3,10), (3,11), (4,7), (4,9), (4,11), (5,8), (5,9), (5,11), (6,9), (6,10), (6,11), (7,11), (8,11), (9,11), (10,11)]
                push!(edges, (sort!([vertices[i][1], vertices[j][1]]), minimum([vertices[i][2], vertices[j][2]])))
            end

            triangles = Vector{Tuple{Vector{Int}, Float64}}([])
            for (i, j, k) in [(1,7,11), (1,8,11), (2,7,11), (2,10,11), (3,8,11), (3,10,11), (4,7,11), (4,9,11), (5,8,11), (5,9,11), (6,9,11), (6,10,11)]
                push!(triangles, ([vertices[i][1], vertices[j][1], vertices[k][1]], minimum([vertices[i][2], vertices[j][2], vertices[k][2]])))
            end

            union!(filtration, Set([([v[1]], v[2]) for v in vertices]))
            union!(filtration, Set(edges))
            union!(filtration, Set(triangles))
        end
    end
    barycenters, sort!(sort!(collect(filtration), by = x -> x[1]), by = x -> length(x[1]))
end