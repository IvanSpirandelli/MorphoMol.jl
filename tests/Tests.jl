module Tests

export run

using Test 
using GeometryBasics

using MorphoMol

include("test_morphometric_approach.jl")
include("test_interface.jl")

function run()
    @testset verbose = true "Morphometric Approach Tests" begin
        run_morphometric_approach_tests()
    end    
    @testset verbose = true "Interface Tests" begin
        run_interface_tests()
    end
end 

end