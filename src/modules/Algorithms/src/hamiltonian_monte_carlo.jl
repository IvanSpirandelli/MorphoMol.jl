struct HamiltonianMonteCarlo{E,EG,IP,DP,EM}
    energy::E
    energy_gradient!::EG
    inner_product::IP
    draw_perturbation!::DP
    exponential_map::EM
    β::Float64              # inverse temperature 
    L::Int                  # number of Leapfrog iterations
    ε::Float64              # Leapfrog stepsize
end

# p = momentum
# x = position
# ∇E = gradient of the energy
# β = inverse temperature
# ε = leapfrog step size
# L = number of leapfrog steps
# grad_energy! = function that computes the gradient of the energy
# exponential_map = exponential map on the lie group
function leapfrog!(p, x, ∇E, β, ε, L, grad_energy!, exponential_map)
    p -= ε*β/2.0 * grad_energy!(∇E, x)
    x = exponential_map(ε * p) * x
    for i in 1:L-1
        p -= ε * β * grad_energy!(∇E, x)
        x = exponential_map(ε * p) * x
    end
    p -= ε * β / 2.0 * grad_energy!(∇E, x)
    x, p
end

function simulate!(hmc::HamiltonianMonteCarlo, x, p, iterations)
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    inner_product, draw_perturbation! = hmc.inner_product, hmc.draw_perturbation!
    exponential_map = hmc.exponential_map
    β, L, ε = hmc.β, hmc.L, hmc.ε

    states = []

    ∇E = zero(p)
    
    x_backup = copy(x)
    E = energy(x)

    accepted_steps = 0
    for _ in 1:iterations
        draw_perturbation!(p)

        H_start = β * E + (1/2)*inner_product(p) # inner_product(p) gives <p, p>
        x, p = leapfrog!(p, x, ∇E, β, ε, L, energy_gradient!, exponential_map)
        E = energy(x)
        E_backup = E

        H_end = β*E + (1/2)*inner_product(p)

        if rand() < min(1, exp(H_start - H_end)) # TODO: drop the min?
            #accept step
            #copyto!(R_backup, R)
            x_backup = copy(x)
            accepted_steps += 1
            push!(states, copy(x))
        else
            # reject step
            E = E_backup
            #copyto!(x, x_backup)
            x = copy(x_backup)
        end
    end
    states, accepted_steps
end