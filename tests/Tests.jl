module Tests

export run

using Test 
using GeometryBasics

using MorphoMol

include("test_morphometric_approach.jl")
include("test_interface.jl")
include("test_energy_calls.jl")

function run()
    @testset verbose = true "Energy Call Tests" begin
        run_energy_call_tests()
    end    
    @testset verbose = true "Morphometric Approach Tests" begin
        run_morphometric_approach_tests()
    end    
    @testset verbose = true "Interface Tests" begin
        run_interface_tests()
    end
end 

end