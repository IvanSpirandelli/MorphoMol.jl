using Rotations, LinearAlgebra

# Your HMC struct is correct
struct HamiltonianMonteCarlo{E,EG,LF}
    energy::E
    energy_gradient!::EG
    leapfrog!::LF
    β::Float64
    L::Int
    ε::Float64
    mass_matrix_diag::Vector{Float64}
end

function muladd!(x, ε, M_inv, p)
    for (i, rb) in enumerate(x)
        # Rotational update using the exponential map
        rot_update_vec = ε .* p[(i-1)*6+1:(i-1)*6+3] .* M_inv[(i-1)*6+1:(i-1)*6+3]
        q_new = rb[1] * exp(Rotations.RotationVec(rot_update_vec...)) # Using RotationVec is slightly more direct

        # Translational update
        trans_update_vec = ε .* p[(i-1)*6+4:i*6] .* M_inv[(i-1)*6+4:i*6]
        t_new = rb[2] .+ trans_update_vec
        
        x[i] = (q_new, t_new)
    end
    # The return is not strictly necessary as x is modified, but it's harmless.
    return x
end

function standard_leapfrog!(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, p::Vector{Float64}, ∇E, β, ε, L, energy_gradient!, M_inv)
    # 1. Half step for momentum
    energy_gradient!(∇E, x)
    p .-= ε * β / 2.0 .* ∇E
    
    # L full steps
    for i in 1:L
        # 2. Full step for position
        muladd!(x, ε, M_inv, p)
        
        # 3. Full step for momentum (except for the last step)
        if i != L
            energy_gradient!(∇E, x)
            p .-= ε * β .* ∇E
        end
    end
    
    # 4. Final half step for momentum
    energy_gradient!(∇E, x)
    p .-= ε * β / 2.0 .* ∇E
    
    # 5. Negate momentum for reversibility
    p .= -p
    
    return x, p 
end

function simulate!(hmc::HamiltonianMonteCarlo, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, iterations::Int)
    # Unpack parameters
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    leapfrog! = hmc.leapfrog!
    β, L, ε = hmc.β, hmc.L, hmc.ε
    mass_matrix_diag = hmc.mass_matrix_diag 
    mass_matrix_inv_diag = 1.0 ./ mass_matrix_diag

    x_current = deepcopy(x)
    
    num_bodies = length(x)
    num_dof = 6 * num_bodies # 3 rotational, 3 translational DOFs per body
    ∇E = zeros(Float64, num_dof)
    
    E_current, _ = energy(x_current)
    
    accepted_steps = 0
    for _ in 1:iterations
        p_current = randn(num_dof) .* sqrt.(mass_matrix_diag)

        # Store current state for potential rejection
        x_proposal = deepcopy(x_current)
        p_proposal = deepcopy(p_current)
        
        # Calculate starting Hamiltonian
        K_current = 0.5 * sum(p_current.^2 .* mass_matrix_inv_diag)
        H_current = β * E_current + K_current

        # Run leapfrog integrator
        leapfrog!(x_proposal, p_proposal, ∇E, β, ε, L, energy_gradient!, mass_matrix_inv_diag)

        # Calculate proposed energy
        E_proposal, _ = energy(x_proposal)
        
        # Calculate ending Hamiltonian
        K_proposal = 0.5 * sum(p_proposal.^2 .* mass_matrix_inv_diag)
        H_proposal = β * E_proposal + K_proposal
        
        # Metropolis-Hastings acceptance step
        acceptance_prob = min(1.0, exp(H_current - H_proposal))
        
        if rand() < acceptance_prob
            x_current = x_proposal
            E_current = E_proposal
            accepted_steps += 1
        else
            # Reject step: Do nothing, the state remains x_current
        end
    end

    return x_current, E_current, accepted_steps / iterations
end

function simulate!(hmc::HamiltonianMonteCarlo, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    # Unpack parameters
    energy, energy_gradient! = hmc.energy, hmc.energy_gradient!
    leapfrog! = hmc.leapfrog!
    β, L, ε = hmc.β, hmc.L, hmc.ε
    mass_matrix_diag = hmc.mass_matrix_diag 
    mass_matrix_inv_diag = 1.0 ./ mass_matrix_diag

    x_current = deepcopy(x)
    
    num_bodies = length(x)
    num_dof = 6 * num_bodies # 3 rotational, 3 translational DOFs per body
    ∇E = zeros(Float64, num_dof)
    
    E_current, measures = energy(x_current)

    total_step_attempts = 1
    add_to_output(merge!(measures, Dict("Es" => E_current, "states" => x, "αs" => total_step_attempts)), output)
    
    while Dates.value(now() - start_time) / 60000.0 < simulation_time_minutes
        total_step_attempts += 1
        p_current = randn(num_dof) .* sqrt.(mass_matrix_diag)

        # Store current state for potential rejection
        x_proposal = deepcopy(x_current)
        p_proposal = deepcopy(p_current)
        
        # Calculate starting Hamiltonian
        K_current = 0.5 * sum(p_current.^2 .* mass_matrix_inv_diag)
        H_current = β * E_current + K_current

        # Run leapfrog integrator
        leapfrog!(x_proposal, p_proposal, ∇E, β, ε, L, energy_gradient!, mass_matrix_inv_diag)

        # Calculate proposed energy
        E_proposal, _ = energy(x_proposal)
        
        # Calculate ending Hamiltonian
        K_proposal = 0.5 * sum(p_proposal.^2 .* mass_matrix_inv_diag)
        H_proposal = β * E_proposal + K_proposal
        
        # Metropolis-Hastings acceptance step
        acceptance_prob = min(1.0, exp(H_current - H_proposal))
        
        if rand() < acceptance_prob
            x_current = x_proposal
            E_current = E_proposal
            add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => total_step_attempts)), output)
        else
            # Reject step: Do nothing, the state remains x_current
        end
    end
    add_to_output(Dict("total_step_attempts" => total_step_attempts), output)
    output
end
