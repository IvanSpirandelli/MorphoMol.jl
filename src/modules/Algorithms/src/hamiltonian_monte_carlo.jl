struct HamiltonianMonteCarlo{E,EG,LF}
    energy::E
    energy_gradient!::EG
    leapfrog!::LF
    β::Float64              # inverse temperature 
    L::Int                  # number of Leapfrog iterations
    ε::Float64              # Leapfrog stepsize
    Σ::Vector{Float64}      # diagonal of covariance matrix
end

# Standard leapfrog algorithm. Assumes that x and p are of the same dimensions.
# p = momentum
# x = position
# ∇E = gradient of the energy
# β = inverse temperature
# ε = leapfrog step size
# L = number of leapfrog steps
# grad_energy! = function that computes the gradient of the energy
function standard_leapfrog!(x, p, ∇E, β, ε, L, energy_gradient!)
    p -= ε * β / 2.0 * energy_gradient!(∇E, x)
    x += ε * p
    for i in 1:L-1
        p -= ε * β * energy_gradient!(∇E, x)
        x += ε * p
    end
    p -= ε * β / 2.0 * energy_gradient!(∇E, x)
    x, p
end

function simulate!(hmc::HamiltonianMonteCarlo, x::Vector{Float64}, iterations::Int)
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    leapfrog! = hmc.leapfrog!
    β, L, ε, Σ = hmc.β, hmc.L, hmc.ε, hmc.Σ

    x_backup = deepcopy(x)
    
    p = zero(x)
    ∇E = zero(x)
    
    E = energy(x)
    
    accepted_steps = 0
    for _ in 1:iterations
        p = randn(length(p)) .* Σ
        Σinv_p = [1.0/e for e in Σ] .* p

        E_backup = E
        H_start = β*E + (1/2) * (p ⋅ Σinv_p)

        x, p = leapfrog!(x, p, ∇E, β, ε, L, energy_gradient!)

        E = energy(x)
        H_end = β*E + (1/2) * (p ⋅ Σinv_p)
        
        if rand() < min(1, exp(H_start - H_end)) 
            #accept step
            copyto!(x_backup, x)
            accepted_steps += 1
        else
            # reject step
            E = E_backup
            copyto!(x, x_backup)
        end
    end
    x, E, accepted_steps/iterations
end


function simulate!(hmc::HamiltonianMonteCarlo, output::SimulationStates, x::Vector{Float64}, iterations::Int)
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    leapfrog! = hmc.leapfrog!
    β, L, ε, Σ = hmc.β, hmc.L, hmc.ε, hmc.Σ

    x_backup = deepcopy(x)
    
    p = zero(x)
    ∇E = zero(x)
    
    E = energy(x)
    add_to_output(x, 0.0, output)

    accepted_steps = 0
    for i in 1:iterations
        p = randn(length(p)) .* Σ
        Σinv_p = [1.0/e for e in Σ] .* p

        E_backup = E
        H_start = β*E + (1/2) * (p ⋅ Σinv_p)

        x, p = leapfrog!(x, p, ∇E, β, ε, L, energy_gradient!)

        E = energy(x)
        H_end = β*E + (1/2) * (p ⋅ Σinv_p)
        
        if rand() < min(1, exp(H_start - H_end)) 
            #accept step
            copyto!(x_backup, x)
            accepted_steps += 1
            add_to_output(x, accepted_steps/iterations, output)
        else
            # reject step
            E = E_backup
            copyto!(x, x_backup)
        end
    end
    output
end


function simulate!(hmc::HamiltonianMonteCarlo, output::MorphometricSimulationOutput, x::Vector{Float64}, iterations::Int)
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    leapfrog! = hmc.leapfrog!
    β, L, ε, Σ = hmc.β, hmc.L, hmc.ε, hmc.Σ

    x_backup = deepcopy(x)
    
    p = zero(x)
    ∇E = zero(x)
    
    E, measures = energy(x)
    add_to_output(x, E, measures, 0.0, output)

    accepted_steps = 0
    
    for i in 1:iterations
        p = randn(length(p)) .* Σ
        Σinv_p = [1.0/e for e in Σ] .* p

        E_backup = E
        H_start = β*E + (1/2) * (p ⋅ Σinv_p)
        x, p = leapfrog!(x, p, ∇E, β, ε, L, energy_gradient!)


        E, measures = energy(x)
        H_end = β*E + (1/2) * (p ⋅ Σinv_p)
        
        if rand() < min(1, exp(H_start - H_end)) 
            #accept step
            copyto!(x_backup, x)
            accepted_steps += 1
            add_to_output(x, E, measures, accepted_steps/i, output)
        else
            # reject step
            E = E_backup
            copyto!(x, x_backup)
        end
    end
    output
end

function simulate!(hmc::HamiltonianMonteCarlo, output::MorphometricSimulationOutput, x::Vector{Float64}, simulation_time_minutes::Float64)
    start_time = now()
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    leapfrog! = hmc.leapfrog!
    β, L, ε, Σ = hmc.β, hmc.L, hmc.ε, hmc.Σ

    x_backup = deepcopy(x)
    
    p = zero(x)
    ∇E = zero(x)
    
    E, measures = energy(x)
    add_to_output(x, E, measures, 0.0, output)

    iterations = 0
    accepted_steps = 0
    while Dates.value(now() - start_time) / 60000.0 < simulation_time_minutes
        iterations += 1
        p = randn(length(p)) .* Σ
        Σinv_p = [1.0/e for e in Σ] .* p

        E_backup = E
        H_start = β*E + (1/2) * (p ⋅ Σinv_p)

        x, p = leapfrog!(x, p, ∇E, β, ε, L, energy_gradient!)

        E, measures = energy(x)
        H_end = β*E + (1/2) * (p ⋅ Σinv_p)
        
        if rand() < min(1, exp(H_start - H_end)) 
            #accept step
            copyto!(x_backup, x)
            accepted_steps += 1
            add_to_output(x, E, measures, accepted_steps/iterations, output)
        else
            # reject step
            E = E_backup
            copyto!(x, x_backup)
        end
    end
    output
end