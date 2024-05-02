module Tests

export run

using Test 

using MorphoMol.Energies

include("test_morphometric_approach.jl")


function run()
    @testset verbose = true "SolSim Tests" begin
        run_morphometric_approach_tests()
    end
end 

end