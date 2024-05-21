struct HamiltonianMonteCarlo{E,EG,IP,DP,LF}
    energy::E
    energy_gradient!::EG
    inner_product::IP
    draw_perturbation::DP
    leapfrog!::LF
    β::Float64              # inverse temperature 
    L::Int                  # number of Leapfrog iterations
    ε::Float64              # Leapfrog stepsize
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

function simulate!(hmc::HamiltonianMonteCarlo, x, iterations)
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    inner_product, draw_perturbation = hmc.inner_product, hmc.draw_perturbation
    leapfrog! = hmc.leapfrog!
    β, L, ε = hmc.β, hmc.L, hmc.ε

    x_backup = deepcopy(x)

    states = [deepcopy(x)]
    
    p = zero(x)
    ∇E = zero(x)
    
    E = energy(x)
    
    accepted_steps = 0
    for _ in 1:iterations
        p = draw_perturbation()

        E_backup = E
        H_start = β*E + (1/2)*inner_product(p)

        x, p = leapfrog!(x, p, ∇E, β, ε, L, energy_gradient!)

        E = energy(x)
        H_end = β*E + (1/2)*inner_product(p)
        
        if rand() < min(1, exp(H_start - H_end)) 
            #accept step
            copyto!(x_backup, x)
            accepted_steps += 1
            push!(states, deepcopy(x))
        else
            # reject step
            E = E_backup
            copyto!(x, x_backup)
        end
    end
    states, accepted_steps
end

struct InputParameters
    algorithm::HamiltonianMonteCarlo
    
end