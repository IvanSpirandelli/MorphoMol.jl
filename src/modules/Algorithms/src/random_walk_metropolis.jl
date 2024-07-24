struct RandomWalkMetropolis{E, P}
    energy::E
    perturbation::P
    β::Float64
end

function simulate!(algorithm::RandomWalkMetropolis, x::Vector{Float64}, iterations::Int)
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E = energy(x)
    accepted_steps = 0
    for _ in 1:iterations
        x_cand = perturbation(x)
        E_backup = energy(x_cand)

        if rand() < exp(-β*(E_backup - E))
            E = E_backup
            accepted_steps += 1
            copy!(x, x_cand)
        end
    end

    return x, accepted_steps/iterations
end

function simulate!(algorithm::RandomWalkMetropolis, output::SimulationStates, x::Vector{Float64}, iterations::Int)
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    x_backup = deepcopy(x)
    E = energy(x)
    add_to_output(x, 0.0, output)

    accepted_steps = 0
    for _ in 1:iterations
        x_backup = perturbation(x)
        E_backup = energy(x_backup)

        if rand() < exp(-β*(E_backup - E))
            E = E_backup
            accepted_steps += 1
            x = deepcopy(x_backup)
            add_to_output(x, accepted_steps/iterations, output)
        end
    end

    return output
end

function simulate!(algorithm::RandomWalkMetropolis, output::MorphometricSimulationOutput, x::Vector{Float64}, simulation_time_minutes::Float64)
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    x_backup = deepcopy(x)

    E, measures = energy(x)
    add_to_output(x, E, measures, 0.0, output)

    accepted_steps = 0
    total_steps = 0
    while Dates.value(now() - start_time) / 60000.0 < simulation_time_minutes
        total_steps += 1
        x_backup = perturbation(x)
        E_backup, measures = energy(x_backup)

        if rand() < exp(-β*(E_backup - E))
            E = E_backup
            accepted_steps += 1
            x = deepcopy(x_backup)
            add_to_output(x, E, measures, accepted_steps/total_steps, output)
        end
    end

    return output
end