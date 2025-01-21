struct ConnectedComponentRandomWalkMetropolis{E, P, ICC}
    energy::E
    perturbation::P
    get_initial_connected_components::ICC
    β::Float64
end

function simulate!(algorithm::ConnectedComponentRandomWalkMetropolis, x::Vector{Float64}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    get_initial_connected_components = algorithm.get_initial_connected_components
    β = algorithm.β

    x_cand = deepcopy(x)
    icc = deepcopy(get_initial_connected_components(x))
    E, measures, _ = energy(icc, 1, x)

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
        i, x_cand = perturbation(x)
        E_cand, measures, ucc = energy(icc, i, x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            accepted_steps += 1
            x = deepcopy(x_cand)
            icc = deepcopy(ucc)
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x, "αs" => accepted_steps/total_steps)), output)
        end
    end

    return output
end

