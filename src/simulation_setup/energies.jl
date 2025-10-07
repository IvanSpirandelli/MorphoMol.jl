using Distributions
using Rotations
using Distances
using Random
using GeometryBasics
using LinearAlgebra

function get_energy(input)
    if input["energy"] == "tasp"
        return (x) -> total_alpha_shape_persistence(x, input["template_centers"], input["persistence_weights"], input["exact_delaunay"])
    elseif input["energy"] == "twasp"
        radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
        return (x) -> total_weighted_alpha_shape_persistence(x, input["template_centers"], radii, input["persistence_weights"], input["exact_delaunay"])
    elseif input["energy"] == "hstasp"
        return (x) -> hs_total_alpha_shape_persistence(x, input["persistence_weights"], input["exact_delaunay"])
    elseif input["energy"] == "fsol"
        radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
        return (x) -> solvation_free_energy_and_measures_in_bounds(x, input["template_centers"], radii, input["rs"], input["prefactors"], input["overlap_jump"], input["overlap_slope"], input["bounds"], input["delaunay_eps"])
    elseif input["energy"] == "cc_fsol"
        return get_connected_component_solvation_free_energy_in_bounds_energy_call(input)
    elseif input["energy"] == "cc_fsol_twasp_interpolated"
        return get_connected_component_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds_energy_call(input, true)
    else    
        return (x) -> 0.0
    end
end

function get_connected_component_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds_energy_call(input, weighted = true)
    mol_type = input["mol_type"]
    rs = input["rs"]
    prefactors = input["prefactors"]
    overlap_jump = input["overlap_jump"]
    overlap_slope = input["overlap_slope"]
    delaunay_eps = input["delaunay_eps"]
    template_centers = input["template_centers"]
    template_radii = input["template_radii"]
    bounds = input["bounds"]
    persistence_weights = input["persistence_weights"]
    exact_delaunay = input["exact_delaunay"]
    μ = input["mu"]

    ssu_energy, ssu_measures = get_single_subunit_energy_and_measures(mol_type, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)

    if input["n_mol"] == 2
        radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
        bol_nmol = (x) -> are_bounding_spheres_overlapping(x, 1, 2, get_bounding_radii(template_centers, template_radii, rs))
        return (x) -> two_mol_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, delaunay_eps, exact_delaunay, ssu_energy, ssu_measures, bol_nmol, μ, weighted)
    else
        bol_nmol = (x, id1, id2) -> are_bounding_spheres_overlapping(x, id1, id2, get_bounding_radii(template_centers, template_radii, rs))
        return (ccs, p_id, x) -> connected_component_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds(ccs, p_id, x, template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, delaunay_eps, exact_delaunay, ssu_energy, ssu_measures, bol_nmol, μ, weighted)
    end
end

function two_mol_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, delaunay_eps, exact_delaunay, ssu_energy, ssu_measures, bol_nmol, μ, compute_weighted::Bool)
    if in_bounds(x, bounds)
        tasp, tasp_measures = compute_weighted ? total_weighted_alpha_shape_persistence(x, template_centers, radii, persistence_weights, exact_delaunay) : total_alpha_shape_persistence(x, template_centers, persistence_weights, exact_delaunay)
        fsol, fsol_measures = solvation_free_energy_and_measures_with_overlap_check(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps, ssu_energy, ssu_measures, bol_nmol)
        μ * fsol + (1 - μ) * tasp, merge!(fsol_measures, tasp_measures)
    else
        return Inf, Dict{String, Any}()
    end
end 

function connected_component_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds(ccs, p_id, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, delaunay_eps, exact_delaunay, ssu_energy, ssu_measures, bol_nmol, μ, compute_weighted::Bool)
    if in_bounds(x, bounds)
        radii = vcat([template_radii for i in 1:length(x)]...)
        tasp, tasp_measures = compute_weighted ? total_weighted_alpha_shape_persistence(x, template_centers, radii, persistence_weights, exact_delaunay) : total_alpha_shape_persistence(x, template_centers, persistence_weights, exact_delaunay)
        fsol, fsol_measures, updated_ccs = connected_component_wise_solvation_free_energy_and_measures(
            ccs,
            p_id,
            x,
            template_centers, 
            template_radii,
            rs, 
            prefactors, 
            overlap_jump,
            overlap_slope,
            delaunay_eps,
            ssu_energy,
            ssu_measures,
            bol_nmol
        )
        μ * fsol + (1 - μ) * tasp, merge!(fsol_measures, tasp_measures), updated_ccs
    else
        return Inf, Dict{String, Any}(), ccs
    end
end 


function get_connected_component_solvation_free_energy_in_bounds_energy_call(input)
    mol_type = input["mol_type"]
    rs = input["rs"]
    prefactors = input["prefactors"]
    overlap_jump = input["overlap_jump"]
    overlap_slope = input["overlap_slope"]
    delaunay_eps = input["delaunay_eps"]
    template_centers = input["template_centers"]
    template_radii = input["template_radii"]
    bounds = input["bounds"]

    ssu_energy, ssu_measures = get_single_subunit_energy_and_measures(mol_type, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)

    if input["n_mol"] == 2
        radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
        bol_nmol = (x) -> are_bounding_spheres_overlapping(x, 1, 2, get_bounding_radii(template_centers, template_radii, rs))
        return (x) -> solvation_free_energy_and_measures_with_overlap_check_in_bounds(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, bounds, delaunay_eps, ssu_energy, ssu_measures, bol_nmol)
    else
        bol_nmol = (x, id1, id2) -> are_bounding_spheres_overlapping(x, id1, id2, get_bounding_radii(template_centers, template_radii, rs))
        return (ccs, p_id, x) -> connected_component_wise_solvation_free_energy_and_measures_in_bounds(ccs, p_id, x, template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, bounds, delaunay_eps, ssu_energy, ssu_measures, bol_nmol)
    end
end

function two_mol_solvation_free_energy_with_total_alpha_shape_persistence_in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, delaunay_eps, exact_delaunay, ssu_energy, ssu_measures, bol_nmol, compute_weighted::Bool)
    if in_bounds(x, bounds)
        tasp, tasp_measures = compute_weighted ? total_weighted_alpha_shape_persistence(x, template_centers, radii, persistence_weights, exact_delaunay) : total_alpha_shape_persistence(x, template_centers, persistence_weights, exact_delaunay)
        fsol, fsol_measures = solvation_free_energy_and_measures_with_overlap_check(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps, ssu_energy, ssu_measures, bol_nmol)
        fsol + tasp, merge!(fsol_measures, tasp_measures)
    else
        return Inf, Dict{String, Any}()
    end
end 

function connected_component_solvation_free_energy_with_total_alpha_shape_persistence_in_bounds(ccs, p_id, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, delaunay_eps, exact_delaunay, ssu_energy, ssu_measures, bol_nmol, compute_weighted::Bool)
    if in_bounds(x, bounds)
        radii = vcat([template_radii for i in 1:length(x)]...)
        tasp, tasp_measures = compute_weighted ? total_weighted_alpha_shape_persistence(x, template_centers, radii, persistence_weights, exact_delaunay) : total_alpha_shape_persistence(x, template_centers, persistence_weights, exact_delaunay)
        fsol, fsol_measures, updated_ccs = connected_component_wise_solvation_free_energy_and_measures(
            ccs,
            p_id,
            x,
            template_centers, 
            template_radii,
            rs, 
            prefactors, 
            overlap_jump,
            overlap_slope,
            delaunay_eps,
            ssu_energy,
            ssu_measures,
            bol_nmol
        )
        fsol + tasp, merge!(fsol_measures, tasp_measures), updated_ccs
    else
        return Inf, Dict{String, Any}(), ccs
    end
end 

function solvation_free_energy(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, delaunay_eps::Float64)
    n_atoms_per_mol = size(template_centers)[2]
    flat_realization = get_flat_realization(x, template_centers)
    Energies.solvation_free_energy(flat_realization, n_atoms_per_mol, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

function get_solvation_free_energy_of_dispersed_configuration(n_mol, mol_type)
    template_centers = MorphoMol.TEMPLATES[mol_type]["template_centers"]
    template_radii = MorphoMol.TEMPLATES[mol_type]["template_radii"]
    pf = MorphoMol.Energies.get_prefactors(1.4, 0.3665)
    n_mol * solvation_free_energy([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], template_centers, template_radii, 1.4, pf, 0.0, 1.1, 100.0)
end

function solvation_free_energy_in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, bounds::Float64, delaunay_eps::Float64)
    if in_bounds(x, bounds)
        n_atoms_per_mol = size(template_centers)[2]
        flat_realization = get_flat_realization(x, template_centers)
        Energies.solvation_free_energy(flat_realization, n_atoms_per_mol, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
    else
        Inf
    end
end

function solvation_free_energy_in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}}, radii::Vector{Vector{Float64}}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, bounds::Float64, delaunay_eps::Float64)
    if in_bounds(x, bounds)
        flat_realization = get_flat_realization(x, template_centers)
        Energies.solvation_free_energy(flat_realization, [size(tc)[2] for tc in template_centers], vcat(radii...), rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
    else
        Inf
    end
end

function solvation_free_energy_and_measures(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}}, radii::Vector{Vector{Float64}}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, delaunay_eps::Float64)
    flat_realization = get_flat_realization(x, template_centers)
    measures = Energies.get_geometric_measures_and_overlap_value(flat_realization, [size(tc)[2] for tc in template_centers], vcat(radii...), rs, overlap_jump, overlap_slope, delaunay_eps)
    sum(measures .* [prefactors; 1.0]), Dict{String,Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5])
end

function solvation_free_energy_and_measures(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, delaunay_eps::Float64)
    n_atoms_per_mol = size(template_centers)[2]
    flat_realization = get_flat_realization(x, template_centers)
    measures = Energies.get_geometric_measures_and_overlap_value(flat_realization, n_atoms_per_mol, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
    sum(measures .* [prefactors; 1.0]), Dict{String,Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5])
end

function in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, bounds::Float64)
    all([all([0.0 <= e <= bounds for e in t]) for (_,t) in x])
end

function solvation_free_energy_and_measures_in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, bounds::Float64, delaunay_eps::Float64)
    if in_bounds(x, bounds)
        n_atoms_per_mol = size(template_centers)[2]
        flat_realization = get_flat_realization(x, template_centers)
        measures = Energies.get_geometric_measures_and_overlap_value(flat_realization, n_atoms_per_mol, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
        sum(measures .* [prefactors; 1.0]), Dict{String,Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5])
    else
        Inf, Dict{String,Any}()
    end
end

function solvation_free_energy_and_measures_in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}}, radii::Vector{Vector{Float64}}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, bounds::Float64, delaunay_eps::Float64)
    if in_bounds(x, bounds)
        flat_realization = get_flat_realization(x, template_centers)
        measures = Energies.get_geometric_measures_and_overlap_value(flat_realization, [size(tc)[2] for tc in template_centers], vcat(radii...), rs, overlap_jump, overlap_slope, delaunay_eps)
        sum(measures .* [prefactors; 1.0]), Dict{String,Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5])
    else
        Inf, Dict{String,Any}()
    end
end

function rotation_and_translation_gradient!(∇E, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, ∇FSol, template_centers)
    n_atoms_per_mol = size(template_centers)[2]
    n_mol = length(x)
    for i in 1:n_mol        
        R = x[i][1]
        ∇E[(i-1) * 6 + 1] = 0.5 * sum([-v[2]*(R[3,:] ⋅ w) + v[3]*(R[2,:] ⋅ w) for (v,w) in [(∇FSol[:,:,i][:,j], template_centers[:,j]) for j in 1:n_atoms_per_mol]])
        ∇E[(i-1) * 6 + 2] = 0.5 * sum([v[1]*(R[3,:] ⋅ w) - v[3]*(R[1,:] ⋅ w) for (v,w) in [(∇FSol[:,:,i][:,j], template_centers[:,j]) for j in 1:n_atoms_per_mol]])
        ∇E[(i-1) * 6 + 3] = 0.5 * sum([-v[1]*(R[2,:] ⋅ w) + v[2]*(R[1,:] ⋅ w) for (v,w) in [(∇FSol[:,:,i][:,j], template_centers[:,j]) for j in 1:n_atoms_per_mol]])
        ∇E[(i-1) * 6 + 4:(i-1) * 6 + 6] = sum([∇FSol[:,j,i] for j in 1:n_atoms_per_mol])
    end
    ∇E
end

function solvation_free_energy_gradient!(∇E, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers, radii, rs, pf, overlap_slope)
    n_atoms_per_mol = size(template_centers)[2]
    n_mol = length(x)
    flat_realization = get_flat_realization(x, template_centers)
    _, dvol, dsurf, dmean, dgauss, dlol = Energies.get_geometric_measures_and_overlap_value_with_derivatives(
        flat_realization,
        n_atoms_per_mol,
        radii,
        rs,
        0.0,
        overlap_slope
    )
    ∇FSol = reshape(pf[1] * dvol + pf[2] * dsurf + pf[3] * dmean + pf[4] * dgauss + dlol, (3, n_atoms_per_mol, n_mol))
    rotation_and_translation_gradient!(∇E, x, ∇FSol, template_centers)
end

function translation_gradient!(∇E, x::Vector{Vector{Float64}}, ∇FSol)
    n_mol = length(x)
    for i in 1:n_mol        
        # Sum the gradients for all atoms in the molecule to get the translational gradient
        ∇E[(i-1) * 3 + 1 : i * 3] = sum(eachcol(∇FSol[:, :, i]))
    end
    ∇E
end

function solvation_free_energy_gradient!(∇E, x::Vector{Vector{Float64}}, radii, rs, pf, overlap_slope)
    n_mol = length(x)
    _, dvol, dsurf, dmean, dgauss, dlol = Energies.get_geometric_measures_and_overlap_value_with_derivatives(
        vcat(x...),
        1,
        radii,
        rs,
        0.0,
        overlap_slope
    )
    ∇FSol = reshape(pf[1] * dvol + pf[2] * dsurf + pf[3] * dmean + pf[4] * dgauss + dlol, (3, 1, n_mol))
    translation_gradient!(∇E, x, ∇FSol,)
end

function hs_total_alpha_shape_persistence(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, persistence_weights::Vector{Float64}, exact_delaunay = false)
    pdgm = Energies.get_alpha_shape_persistence_diagram(collect(eachcol(reshape(x, (3,length(x)÷3)))), exact_delaunay)
    p0 = Energies.get_total_persistence(pdgm[1])
    p1 = Energies.get_total_persistence(pdgm[2])
    p2 = Energies.get_total_persistence(pdgm[3])
    λ0, λ1, λ2 = persistence_weights
    λ0 * p0 + λ1 * p1 + λ2 * p2, Dict{String, Any}("P0s" => p0, "P1s" => p1, "P2s" => p2)
end

function total_alpha_shape_persistence(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64}, persistence_weights::Vector{Float64}, exact_delaunay = false)
    pdgm = Energies.get_alpha_shape_persistence_diagram(get_point_vector_realization(x, template_centers), exact_delaunay)
    p0 = Energies.get_total_persistence(pdgm[1])
    p1 = Energies.get_total_persistence(pdgm[2])
    p2 = Energies.get_total_persistence(pdgm[3])
    λ0, λ1, λ2 = persistence_weights
    λ0 * p0 + λ1 * p1 + λ2 * p2, Dict{String, Any}("P0s" => p0, "P1s" => p1, "P2s" => p2)
end

function total_weighted_alpha_shape_persistence(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64}, radii::Vector{Float64}, persistence_weights::Vector{Float64}, exact_delaunay = false)
    pdgm = Energies.get_weighted_alpha_shape_persistence_diagram(get_point_vector_realization(x, template_centers), radii, exact_delaunay)
    p0 = Energies.get_total_persistence(pdgm[1])
    p1 = Energies.get_total_persistence(pdgm[2])
    p2 = Energies.get_total_persistence(pdgm[3])
    λ0, λ1, λ2 = persistence_weights
    λ0 * p0 + λ1 * p1 + λ2 * p2, Dict{String, Any}("P0s" => p0, "P1s" => p1, "P2s" => p2)
end

function solvation_free_energy_with_persistence(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, persistence_weights::Vector{Float64}, delaunay_eps::Float64)
    n_atoms_per_mol = size(template_centers)[2]
    flat_realization = get_flat_realization(x, template_centers)
    points = Vector{Vector{Float64}}([e for e in eachcol(reshape(flat_realization, (3, Int(length(flat_realization) / 3))))])
    pdgm = Energies.get_alpha_shape_persistence_diagram(points)
    p0 = Energies.get_total_persistence(pdgm[1])
    p1 = Energies.get_total_persistence(pdgm[2])
    p2 = Energies.get_total_persistence(pdgm[3])
    λ0, λ1, λ2 = persistence_weights
    measures = Energies.get_geometric_measures_and_overlap_value(flat_realization, n_atoms_per_mol, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
    sum(measures .* [prefactors; [1.0]]) + λ0 * p0 + λ1 * p1 + λ2 * p2, Dict{String, Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5], "P0s" => p0, "P1s" => p1, "P2s" => p2)
end
