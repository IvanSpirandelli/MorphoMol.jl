struct MixedEnergyRandomWalkMetropolis{EA, EB, P}
    energy_a::EA
    energy_b::EB
    perturbation::P
    β_a::Float64
    β_b::Float64
end

function simulate!(algorithm::MixedEnergyRandomWalkMetropolis, x::Vector{Float64}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy_a = algorithm.energy_a
    energy_b = algorithm.energy_b
    β_a = algorithm.β_a
    β_b = algorithm.β_b
    perturbation = algorithm.perturbation

    E_a, measures_a = energy_a(x)
    E_b, measures_b = energy_b(x)
    add_to_output(merge!(measures_a, Dict("Es_a" => E_a, "states" => x, "αs_a" => 0.0)), output)
    add_to_output(merge!(measures_b, Dict("Es_b" => E_b, "states" => x, "αs_b" => 0.0)), output)

    x_backup = deepcopy(x)
    accepted_steps_a = 0
    accepted_steps_b = 0
    total_steps_a = 1
    total_steps_b = 1
    while Dates.value(now() - start_time) / 60000.0 < simulation_time_minutes
        total_steps = total_steps_a + total_steps_b
        if total_steps % 2 == 0
            x_backup = perturbation(x)
            E_backup_a, measures_a = energy_a(x_backup)

            if rand() < exp(-β_a*(E_backup_a - E_a))
                E_a = E_backup_a
                accepted_steps_a += 1
                x = deepcopy(x_backup)
                add_to_output(merge!(measures_a, Dict("Es_a" => E_a, "states" => x, "αs_a" => accepted_steps_a/total_steps_a)), output)
                E_b, measures_b = energy_b(x)
                add_to_output(merge!(measures_b, Dict("Es_b" => E_b)), output)
            end
            total_steps_a += 1
        elseif total_steps % 2 == 1
            x_backup = perturbation(x)
            E_backup_b, measures_b = energy_b(x_backup)

            if rand() < exp(-β_b*(E_backup_b - E_b))
                E_b = E_backup_b
                accepted_steps_b += 1
                x = deepcopy(x_backup)
                add_to_output(merge!(measures_b,Dict("Es_b" => E_b, "states" => x, "αs_b" => accepted_steps_b/total_steps_b)), output)
                E_a, measures_a = energy_a(x)
                add_to_output(merge!(measures_a, Dict("Es_a" => E_a)), output)
            end
            total_steps_b += 1
        end
    end

    return output
end
