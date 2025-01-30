function get_perturbation(input)
    σ_r = input["σ_r"]
    σ_t = input["σ_t"]
    n_mol = input["n_mol"]
    if input["perturbation"] == "single_random"
        return (x) -> perturb_single_randomly_chosen(x, σ_r, σ_t)
    elseif input["perturbation"] == "single_random_get_index"
        return (x) -> get_index_and_perturb_single_randomly_chosen(x, σ_r, σ_t)
    elseif input["perturbation"] == "all"
        Σ = vcat([[σ_r, σ_r, σ_r, σ_t, σ_t, σ_t] for _ in 1:n_mol]...)
        return (x) -> perturb_all(x, Σ)
    end
end

function perturb_all(x, Σ) 
    x .+ (randn(length(x)) .* Σ)
end

function perturb_single_randomly_chosen(x, σ_r, σ_t)
    x_cand = deepcopy(x)
    i  = rand(0:(length(x)÷6)-1)
    x_cand[(i*6)+1:(i*6)+6] = x_cand[(i*6)+1:(i*6)+6] .+ (randn(6) .* [σ_r, σ_r, σ_r, σ_t, σ_t, σ_t])
    x_cand
end

function get_index_and_perturb_single_randomly_chosen(x, σ_r, σ_t)
    x_cand = deepcopy(x)
    i  = rand(0:(length(x)÷6)-1)
    x_cand[(i*6)+1:(i*6)+6] = x_cand[(i*6)+1:(i*6)+6] .+ (randn(6) .* [σ_r, σ_r, σ_r, σ_t, σ_t, σ_t])
    i+1, x_cand
end