function get_initialization(input)
    n_mol = input["n_mol"]
    bounds = input["bounds"]
    if input["initialization"] == "random"
        return (x) -> get_initial_state(n_mol, bounds)
    elseif input["initialization"] == "random_only_translations"
        return (x) -> get_initial_state_only_translations(n_mol, bounds)
    end
end

function get_initial_state(n_mol::Int, bounds::Float64)
    vcat([
        [rand(Uniform(0.0, 2*pi)), rand(Uniform(0.0, 2*pi)), rand(Uniform(0.0, 2*pi)), 
        rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds))] 
        for i in 1:n_mol]...);
end

function get_initial_state_only_translations(n_mol::Int, bounds::Float64)
    vcat([
        [rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds))] 
        for i in 1:n_mol]...);
end