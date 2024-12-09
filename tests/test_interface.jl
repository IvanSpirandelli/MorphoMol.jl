function run_interface_tests()
    @testset verbose = true "Interface" begin
        @testset verbose = true "Simplex Count" begin
            test_simplex_count()
        end
    end
end

function test_simplex_count()
    points = [Point3f([0.0, 0.0, 0.0]), Point3f([2.0, 0.0, 0.0]), Point3f([1.0, 0.0, 2.4]), Point3f([1.0, 3.0, 1.2]), Point3f([4.0, 2.0, 1.8])]
    n_atoms_per_mol = 3
    mc_tets = get_multichromatic_tetrahedra(points, n_atoms_per_mol)
    bcs, filtration = get_barycentric_subdivision_and_filtration(points, mc_tets, n_atoms_per_mol)
    @test length([e for e in filtration if length(e[1]) == 1]) == 13
    @test length([e for e in filtration if length(e[1]) == 2]) == 26
    @test length([e for e in filtration if length(e[1]) == 3]) == 14

    points = [Point3f([0.1, 0.1, 0.5]), Point3f([2.0, 0.0, 0.2]), Point3f([1.0, 2.5, 0.0]), Point3f([-2.0, 2.0, 0.4]), Point3f([1.0, 1.4, 2.0])]
    n_atoms_per_mol = 2
    mc_tets = Matrix{Int64}([1 2 3 5; 1 3 4 5])
    bcs, filtration = get_barycentric_subdivision_and_filtration(points, mc_tets, n_atoms_per_mol)
    @test length([e for e in filtration if length(e[1]) == 1]) == 16 #twice the number of a single tricolored tetrahedron
    @test length([e for e in filtration if length(e[1]) == 2]) == 35 #twice the number of a single tricolored tetrahedron - 3 shared edges
    @test length([e for e in filtration if length(e[1]) == 3]) == 20 #twice the number of a single tricolored tetrahedron
end