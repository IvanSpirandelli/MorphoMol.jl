struct RandomWalkMetropolis{E, P}
    energy::E
    perturbation::P
    β::Float64
end

function simulate!(algorithm::RandomWalkMetropolis, x_init::Vector{Float64}, iterations::Int)
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β
    x = deepcopy(x_init)
    x_cand = deepcopy(x)

    E = energy(x)
    accepted_steps = 0
    for _ in 1:iterations
        x_cand = perturbation(x)
        E_cand = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            accepted_steps += 1
            x = deepcopy(x_cand)
        end
    end

    return x, E, accepted_steps/iterations
end

function simulate!(algorithm::RandomWalkMetropolis, x::Vector{Float64}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    x_cand = deepcopy(x)

    E, measures = energy(x)

    accepted_steps = 0
    total_steps = 0
    # If this is true it means the simulation is a continuation of a previous one
    if haskey(output, "αs") && length(output["αs"]) > 1 
        accepted_steps = length(output["Es"])
        total_steps = Int(round(length(output["Es"]) / output["αs"][end]))
    elseif haskey(output, "αs") && length(output["αs"]) == 1 # This means we reset to acceptance rate to see the one of the last continuation.
        nothing
    else
        add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => 0.0)), output)
    end

    while Dates.value(now() - start_time) / 60000.0 < simulation_time_minutes
        total_steps += 1
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            accepted_steps += 1
            x = deepcopy(x_cand)
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x, "αs" => accepted_steps/total_steps)), output)
        end
    end

    return output
end

