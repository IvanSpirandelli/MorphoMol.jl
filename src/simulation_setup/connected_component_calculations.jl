using Graphs

function get_single_subunit_energy_and_measures(mol_type, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
    template_centers = TEMPLATES[mol_type]["template_centers"]
    template_radii = TEMPLATES[mol_type]["template_radii"]
    solvation_free_energy_and_measures(fill(0.0, 6), template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

function get_bounding_radius(mol_type)
    if mol_type == "6r7m"
        return 43.0
    else
        @assert false
    end
end

function are_bounding_spheres_overlapping(x::Vector{Float64}, id_one::Int, id_two::Int, bounding_radius::Float64)
    t_one = [x[(id_one-1)*6+4], x[(id_one-1)*6+5], x[(id_one-1)*6+6]]
    t_two = [x[(id_two-1)*6+4], x[(id_two-1)*6+5], x[(id_two-1)*6+6]]
    euclidean(t_one, t_two) < 2.0 * bounding_radius
end

#Use with standard RWM and 2 subunits only
function solvation_free_energy_and_measures_with_overlap_check_in_bounds(
    x::Vector{Float64},
    template_centers::Matrix{Float64}, 
    radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    delaunay_eps::Float64,
    single_subunit_energy::Float64,
    single_subunit_measures::Dict{String, Any},
    molecule_boundary_overlap_check::Function
    )
    if in_bounds(x, bounds)
        solvation_free_energy_and_measures_with_overlap_check(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps, single_subunit_energy, single_subunit_measures, molecule_boundary_overlap_check)
    else
        Inf, Dict{String, Any}()
    end
end

function solvation_free_energy_and_measures_with_overlap_check(
    x::Vector{Float64},
    template_centers::Matrix{Float64}, 
    radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    single_subunit_energy::Float64,
    single_subunit_measures::Dict{String, Any},
    molecule_boundary_overlap_check::Function
    )
    if !molecule_boundary_overlap_check(x)
        return 2*single_subunit_energy, Dict{String, Any}(k => v.*2 for (k,v) in single_subunit_measures)
    end
    solvation_free_energy_and_measures(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

#Stuff for 3 or more subunits
get_sub_state(x, indices) = vcat([x[(i-1)*6+1:(i-1)*6+6] for i in indices]...)

function construct_overlap_graph(x::Vector{Float64}, molecule_boundary_overlap_check)
    n = Int(length(x) / 6)
    graph = SimpleGraph(n)
    for i in 1:n
        for j = i+1:n
            if molecule_boundary_overlap_check(x, i, j)
                add_edge!(graph, i, j)
            end
        end
    end
    graph
end

function connected_component_wise_solvation_free_energy_and_measures_in_bounds(last_iteration_ccs_energies_and_measures::Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}},
    transformed_index::Int,
    x::Vector{Float64},
    template_centers::Matrix{Float64}, 
    template_radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    delaunay_eps::Float64,
    single_subunit_energy::Float64,
    single_subunit_measures::Dict{String, Any},
    molecule_boundary_overlap_check::Function,
    )
    if in_bounds(x, bounds)
        connected_component_wise_solvation_free_energy_and_measures(last_iteration_ccs_energies_and_measures, transformed_index, x, template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps, single_subunit_energy, single_subunit_measures, molecule_boundary_overlap_check)
    else
        Inf, Dict{String, Any}(), Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}}()
    end
end

function connected_component_wise_solvation_free_energy_and_measures(
    last_iteration_ccs_energies_and_measures::Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}},
    transformed_index::Int,
    x::Vector{Float64},
    template_centers::Matrix{Float64}, 
    template_radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    single_subunit_energy::Float64,
    single_subunit_measures::Dict{String, Any},
    molecule_boundary_overlap_check::Function,
    )

    graph = construct_overlap_graph(x, molecule_boundary_overlap_check)
    ccs = connected_components(graph)
    indexed_cc = [e for e in connected_components(graph) if transformed_index in e][1]
    ccs_energies_and_measures = Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}}()

    for cc in ccs
        if length(cc) == 1
            ccs_energies_and_measures[cc] = single_subunit_energy, single_subunit_measures
            #println("Setting to single sub unit energy and measure")
            # Setting to single submeasure
        elseif !(cc in keys(last_iteration_ccs_energies_and_measures)) || cc == indexed_cc
            ccs_energies_and_measures[cc] = solvation_free_energy_and_measures(get_sub_state(x, cc), template_centers, vcat([template_radii for _ in 1:length(cc)]...), rs, prefactors, overlap_jump, overlap_slope, delaunay_eps) 
            #println("Recalculating Changed!", ccs_energies_and_measures[cc])
        else
            ccs_energies_and_measures[cc] = last_iteration_ccs_energies_and_measures[cc]
            #println("Setting to Previous!", ccs_energies_and_measures[cc])
        end
    end
    f_sol = 0.0
    measures = Dict{String, Any}(k => 0.0 for k in keys(first(ccs_energies_and_measures)[2][2]))
    for comb_val in values(ccs_energies_and_measures)
        f_sol += comb_val[1]
        for (k,v) in comb_val[2]
            measures[k] += v
        end
    end
    f_sol, measures, ccs_energies_and_measures
end

function get_initial_connected_component_energies(
    x::Vector{Float64},
    template_centers::Matrix{Float64}, 
    template_radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    molecule_boundary_overlap_check::Function
    )

    graph = construct_overlap_graph(x, molecule_boundary_overlap_check)
    ccs = connected_components(graph)
    ccs_energies_and_measures = Dict{Vector{Int}, Tuple{Float64, Dict{String, Any}}}()

    for cc in ccs
        ccs_energies_and_measures[cc] = solvation_free_energy_and_measures(get_sub_state(x, cc), template_centers, vcat([template_radii for _ in 1:length(cc)]...), rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)    
    end
    ccs_energies_and_measures
end