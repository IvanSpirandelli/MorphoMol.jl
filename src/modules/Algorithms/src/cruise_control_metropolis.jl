struct CruiseControlMetropolis{ME, CE, P}
    main_energy::ME
    cruise_energy::CE
    perturbation::P
    β::Float64
    cruise_trail::Int
    cruise_impact_target::Float64
    minimum_λ::Float64
end

function get_average_energy_difference(Es)
    if length(Es) < 2
        return 0.0
    end
    sum([Es[i] - Es[i-1] for i in 2:length(Es)])/length(Es)
end

function simulate!(algorithm::CruiseControlMetropolis, x::Vector{Float64}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    main_energy = algorithm.main_energy
    cruise_energy = algorithm.cruise_energy
    perturbation = algorithm.perturbation
    β = algorithm.β
    cruise_trail = algorithm.cruise_trail
    cruise_impact_target = algorithm.cruise_impact_target
    minimum_λ = algorithm.minimum_λ

    x_cand = deepcopy(x)

    main_E, main_measures = main_energy(x)
    cruise_E, cruise_measures = cruise_energy(x)

    λ = minimum_λ
    E = main_E + λ * cruise_E

    accepted_steps = 0
    total_steps = 0

    add_to_output(merge!(main_measures, Dict("main_Es" => main_E, "states" => x, "αs" => 0.0)), output)
    add_to_output(merge!(cruise_measures, Dict("cruise_Es" => cruise_E, "λs" => λ)), output)

    while Dates.value(now() - start_time) / 60000.0 < simulation_time_minutes
        total_steps += 1

        x_cand = perturbation(x)

        main_diff = get_average_energy_difference(output["main_Es"][maximum([1, length(output["main_Es"]) - cruise_trail]):end])
        cruise_diff = get_average_energy_difference(output["cruise_Es"][maximum([1, length(output["cruise_Es"]) - cruise_trail]):end])

        if total_steps < 1000
            λ = minimum_λ
        else
            λ = cruise_diff == 0.0 ? minimum_λ : maximum([minimum_λ, main_diff * cruise_impact_target / cruise_diff])
        end

        main_E_cand, main_measures = main_energy(x_cand)
        cruise_E_cand, cruise_measures = cruise_energy(x_cand)

        E_cand = main_E_cand + λ * cruise_E_cand

        if rand() < exp(-β*(E_cand - E)) 
            accepted_steps += 1
            E = E_cand
            x = deepcopy(x_cand)
            add_to_output(merge!(main_measures, Dict("main_Es" => main_E_cand, "states" => x, "αs" => accepted_steps/total_steps)), output)
            add_to_output(merge!(cruise_measures, Dict("cruise_Es" => cruise_E_cand, "λs" => λ)), output)
        end
    end
    add_to_output(Dict("αs" => accepted_steps/total_steps), output)
    println("accepted_steps: ", accepted_steps)
    println("total_steps: ", total_steps)
    return output
end