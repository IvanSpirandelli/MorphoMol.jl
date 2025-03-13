using Distances

function get_initialization(input, non_overlapping = false)
    n_mol = input["n_mol"]
    bounds = input["bounds"]
    if input["initialization"] == "random"
        if non_overlapping  
            return () -> get_non_overlapping_initial_state(input)
        else
            return () -> get_initial_state(n_mol, bounds)
        end
    elseif input["initialization"] == "random_only_translations"
        return () -> get_initial_state_only_translations(n_mol, bounds)
    end
end

function get_non_overlapping_initial_state(input, max_count = 25)
    n_mol = input["n_mol"]
    bounds = input["bounds"]
    x = get_initial_state(n_mol, bounds)
    attempt = 0
    while is_overlapping_state(x, input)
        x = get_initial_state(n_mol, bounds)
        attempt += 1
        if attempt >= max_count
            println("Maximum attempts to find non overlapping initial state reached.")
            @assert false 
        end
    end
    return x
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

function is_overlapping_state(x, input)
    realz = get_flat_realization(x, input["template_centers"])
    measures = get_geometric_measures_and_overlap_value(
        realz, 
        length(input["template_radii"]), 
        vcat([input["template_radii"] for i in 1:input["n_mol"]]), 
        0.0, 
        1.0,
        0.0, 
        100.0
    )
    return !isapprox(0.0, measures[5])
end