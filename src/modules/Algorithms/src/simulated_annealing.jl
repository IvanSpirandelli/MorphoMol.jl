struct SimulatedAnnealing{E, P, TR}
    energy::E
    perturbation::P
    temperature_reduction::TR
end

function simulate!(algorithm::SimulatedAnnealing, x_init::Vector{Float64}, iterations::Int, output::Dict)
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    temperature_reduction = algorithm.temperature_reduction
    
    T = temperature_reduction(0)
    x = deepcopy(x_init)
    x_cand = deepcopy(x)
    E = energy(x)
    add_to_output(Dict("Es" => E, "states" => x), output)
    accepted_steps = 0
    for i in 1:iterations
        x_cand = perturbation(x)
        E_backup = energy(x_cand)
        if rand() < exp(-(1.0/T)*(E_backup - E))
            E = E_backup
            accepted_steps += 1
            x = deepcopy(x_cand)
            add_to_output(Dict("Es" => E, "states" => x), output)
        end

        T = temperature_reduction(i)
    end

    return x, E, accepted_steps/iterations
end