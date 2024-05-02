struct HamiltonianMonteCarlo{E,EG,IP,DP}
    energy::E
    energy_gradient!::EG
    inner_product::IP
    draw_perturbation!::DP
    β::Float64              # inverse temperature 
    L::Int                  # number of Leapfrog iterations
    ε::Float64              # Leapfrog stepsize
end

# p = momentum
# R = rotation (current position)
# ∇E = gradient of the energy
# β = inverse temperature
# ε = leapfrog step size
# L = number of leapfrog steps
# grad_energy! = function that computes the gradient of the energy
function leapfrog!(p, R, ∇E, β, ε, L, grad_energy!)
    p -= ε*β/2.0 * grad_energy!(∇E, R)
    R = exp(Rotations.RotationVecGenerator((ε * p)...)) * R
    for i in 1:L-1
        p -= ε * β * grad_energy!(∇E, R)
        R = exp(Rotations.RotationVecGenerator((ε * p)...)) * R
    end
    p -= ε * β / 2.0 * grad_energy!(∇E, R)
    R, p
end

function simulate!(hmc::HamiltonianMonteCarlo, R, p, iterations)
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    inner_product, draw_perturbation! = hmc.inner_product, hmc.draw_perturbation!
    β, L, ε = hmc.β, hmc.L, hmc.ε

    states = []

    ∇E = zero(p)
    
    R_backup = copy(R)
    E = energy(R)

    accepted_steps = 0
    for _ in 1:iterations
        draw_perturbation!(p)

        H_start = β * E + (1/2)*inner_product(p) # inner_product(p) gives <p, p>
        R, p = leapfrog!(p, R, ∇E, β, ε, L, energy_gradient!)
        E = energy(R)
        E_backup = E

        H_end = β*E + (1/2)*inner_product(p)

        if rand() < min(1, exp(H_start - H_end)) # TODO: drop the min?
            #accept step
            #copyto!(R_backup, R)
            R_backup = copy(R)
            accepted_steps += 1
            push!(states, copy(R))
        else
            # reject step
            E = E_backup
            #copyto!(R, R_backup)
            R = copy(R_backup)
        end
    end
    states, accepted_steps
end