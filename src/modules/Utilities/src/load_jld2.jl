function load_jld2(file)
    mutable struct SimulationOutput
        states::Vector{Vector{Float64}}
        Es::Vector{Float64}
        measures::Vector{Vector{Float64}}
        αs::Vector{Float32}
    end

    mutable struct SimulationStates
        states::Vector{Vector{Float64}}
        αs::Vector{Float32}
    end

    mutable struct MorphometricSimulationOutput
        states::Vector{Vector{Float64}}
        Es::Vector{Float64}
        Vs::Vector{Float32}
        As::Vector{Float32}
        Ms::Vector{Float32}
        Xs::Vector{Float32}
        OLs::Vector{Float32}
        αs::Vector{Float32}
    end
end