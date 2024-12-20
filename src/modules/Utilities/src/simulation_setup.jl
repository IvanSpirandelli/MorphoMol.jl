function get_initial_state(n_mol::Int, bounds::Float64)
    vcat([
        [rand(Uniform(0.0, 2*pi)), rand(Uniform(0.0, 2*pi)), rand(Uniform(0.0, 2*pi)), 
        rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds))] 
        for i in 1:n_mol]...);
end

perturb_all(x, Σ) = x .+ (randn(length(x)) .* Σ)

function perturb_single_randomly_chosen(x, σ_r, σ_t)
    x_cand = deepcopy(x)
    i  = rand(0:(length(x)÷6)-1)
    x_cand[(i*6)+1:(i*6)+6] = x_cand[(i*6)+1:(i*6)+6] .+ (randn(6) .* [σ_r, σ_r, σ_r, σ_t, σ_t, σ_t])
    x_cand
end

function solvation_free_energy(x::Vector{Float64}, template_mol::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, delaunay_eps::Float64)
    n_atoms_per_mol = size(template_mol)[2]
    flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_mol)
    MorphoMol.Energies.solvation_free_energy(flat_realization, n_atoms_per_mol, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

function solvation_free_energy_and_measures(x::Vector{Float64}, template_mol::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, delaunay_eps::Float64)
    n_atoms_per_mol = size(template_mol)[2]
    flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_mol)
    measures = MorphoMol.Energies.get_geometric_measures_and_overlap_value(flat_realization, n_atoms_per_mol, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
    sum(measures .* [prefactors; 1.0]), Dict{String,Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5])
end

function in_bounds(x::Vector{Float64}, bounds::Float64)
    all(0.0 <= e && e <= bounds for e in x[4:6:end]) && all(0.0 <= e && e <= bounds for e in x[5:6:end]) && all(0.0 <= e && e <= bounds for e in x[6:6:end])
end

function solvation_free_energy_and_measures_in_bounds(x::Vector{Float64}, template_mol::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, bounds::Float64, delaunay_eps::Float64)
    if in_bounds(x, bounds)
        n_atoms_per_mol = size(template_mol)[2]
        flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_mol)
        measures = MorphoMol.Energies.get_geometric_measures_and_overlap_value(flat_realization, n_atoms_per_mol, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
        sum(measures .* [prefactors; 1.0]), Dict{String,Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5])
    else
        Inf, Dict{String,Any}()
    end
end

function rotation_and_translation_gradient!(∇E, x, ∇FSol, template_mol)
    n_atoms_per_mol = size(template_mol)[2]
    n_mol = length(x) ÷ 6
    for i in 1:n_mol        
        R = exp(Rotations.RotationVecGenerator(x[(i-1)*6 + 1:(i-1)*6 + 3]...))
        ∇E[(i-1) * 6 + 1] = 0.5 * sum([-v[2]*(R[3,:] ⋅ w) + v[3]*(R[2,:] ⋅ w) for (v,w) in [(∇FSol[:,:,i][:,j], template_mol[:,j]) for j in 1:n_atoms_per_mol]])
        ∇E[(i-1) * 6 + 2] = 0.5 * sum([v[1]*(R[3,:] ⋅ w) - v[3]*(R[1,:] ⋅ w) for (v,w) in [(∇FSol[:,:,i][:,j], template_mol[:,j]) for j in 1:n_atoms_per_mol]])
        ∇E[(i-1) * 6 + 3] = 0.5 * sum([-v[1]*(R[2,:] ⋅ w) + v[2]*(R[1,:] ⋅ w) for (v,w) in [(∇FSol[:,:,i][:,j], template_mol[:,j]) for j in 1:n_atoms_per_mol]])
        ∇E[(i-1) * 6 + 4:(i-1) * 6 + 6] = sum([∇FSol[:,j,i] for j in 1:n_atoms_per_mol])
    end
    ∇E
end

function solvation_free_energy_gradient!(∇E, x, template_mol, radii, rs, pf, overlap_slope)
    n_atoms_per_mol = size(template_mol)[2]
    n_mol = length(x) ÷ 6
    flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_mol)
    _, dvol, dsurf, dmean, dgauss, dlol = MorphoMol.Energies.get_geometric_measures_and_overlap_value_with_derivatives(
        flat_realization,
        n_atoms_per_mol,
        radii,
        rs,
        0.0,
        overlap_slope
    )
    ∇FSol = reshape(pf[1] * dvol + pf[2] * dsurf + pf[3] * dmean + pf[4] * dgauss + dlol, (3, n_atoms_per_mol, n_mol))
    rotation_and_translation_gradient!(∇E, x, ∇FSol, template_mol)
end

function persistence(x::Vector{Float64}, template_centers::Matrix{Float64}, persistence_weights::Vector{Float64})
    flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_centers)
    points = Vector{Vector{Float64}}([e for e in eachcol(reshape(flat_realization, (3, Int(length(flat_realization) / 3))))])
    pdgm = MorphoMol.Energies.get_alpha_shape_persistence_diagram(points)
    p0 = MorphoMol.Energies.get_total_persistence(pdgm[1], persistence_weights[1])
    p1 = MorphoMol.Energies.get_total_persistence(pdgm[2], persistence_weights[2])
    p2 = MorphoMol.Energies.get_total_persistence(pdgm[3], persistence_weights[3])
    p0 + p1 + p2, Dict{String, Any}("P0s" => p0, "P1s" => p1, "P2s" => p2)
end

function persistence_with_diagram(x::Vector{Float64}, template_centers::Matrix{Float64}, persistence_weights::Vector{Float64})
    flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_centers)
    points = Vector{Vector{Float64}}([e for e in eachcol(reshape(flat_realization, (3, Int(length(flat_realization) / 3))))])
    pdgm = MorphoMol.Energies.get_alpha_shape_persistence_diagram(points)
    pdgm = [pdgm[1], pdgm[2], pdgm[3]]
    p0 = MorphoMol.Energies.get_total_persistence(pdgm[1], persistence_weights[1])
    p1 = MorphoMol.Energies.get_total_persistence(pdgm[2], persistence_weights[2])
    p2 = MorphoMol.Energies.get_total_persistence(pdgm[3], persistence_weights[3])
    p0 + p1 + p2, Dict{String, Any}("P0" => p0, "P1" => p1, "P2" => p2, "PDGMs"  => pdgm)
end

function interface_persistence(x::Vector{Float64}, template_centers::Matrix{Float64}, persistence_weights::Vector{Float64})
    n_atoms_per_mol = size(template_centers)[2]
    points = MorphoMol.Utilities.get_point_vector_realization(x, template_centers)
    idgm = MorphoMol.Energies.get_interface_persistence_diagram(points, n_atoms_per_mol)
    p0 = MorphoMol.Energies.get_death_by_birth_persistence(idgm[1], persistence_weights[1])
    p1 = MorphoMol.Energies.get_death_by_birth_persistence(idgm[2], persistence_weights[2])
    p0 + p1, Dict{String, Any}("P0s" => p0, "P1s" => p1)
end

function solvation_free_energy_with_persistence(x::Vector{Float64}, template_centers::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, persistence_weights::Vector{Float64}, delaunay_eps::Float64)
    n_atoms_per_mol = size(template_centers)[2]
    flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_centers)
    points = Vector{Vector{Float64}}([e for e in eachcol(reshape(flat_realization, (3, Int(length(flat_realization) / 3))))])
    pdgm = MorphoMol.Energies.get_alpha_shape_persistence_diagram(points)
    pdgm = [pdgm[1], pdgm[2], pdgm[3]]
    p0 = MorphoMol.Energies.get_total_persistence(pdgm[1], persistence_weights[1])
    p1 = MorphoMol.Energies.get_total_persistence(pdgm[2], persistence_weights[2])
    p2 = MorphoMol.Energies.get_total_persistence(pdgm[3], persistence_weights[3])
    measures = MorphoMol.Energies.get_geometric_measures_and_overlap_value(flat_realization, n_atoms_per_mol, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
    sum(measures .* [prefactors; [1.0]]) + p0 + p1 + p2, Dict{String, Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5], "P0s" => p0, "P1s" => p1, "P2s" => p2)
end


function solvation_free_energy_with_interface_persistence_and_measures_in_bounds(x::Vector{Float64}, template_centers::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, persistence_weights::Vector{Float64}, bounds::Float64, delaunay_eps::Float64)
    if in_bounds(x, bounds)
        n_atoms_per_mol = size(template_centers)[2]
        flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_centers)
        points = Vector{Vector{Float64}}([e for e in eachcol(reshape(flat_realization, (3, Int(length(flat_realization) / 3))))])
        idgm = MorphoMol.Energies.get_interface_persistence_diagram(points, n_atoms_per_mol)
        p0 = MorphoMol.Energies.get_death_by_birth_persistence(idgm[1], persistence_weights[1])
        p1 = MorphoMol.Energies.get_death_by_birth_persistence(idgm[2], persistence_weights[2])
        measures = MorphoMol.Energies.get_geometric_measures_and_overlap_value(flat_realization, n_atoms_per_mol, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
        sum(measures .* [prefactors; [1.0]]) + p0 + p1, Dict{String, Any}("Vs" => measures[1], "As" => measures[2], "Cs" => measures[3], "Xs" => measures[4], "OLs" => measures[5], "P0s" => p0, "P1s" => p1)
    else
        Inf, Dict{String, Any}()
    end
end

function calculate_T0(Es, target_acceptance_rate)
    transitions = []
    for i in 1:length(Es)-1
        if Es[i] > Es[i+1]
            push!(transitions, Es[i])
            push!(transitions, Es[i+1])
        end
    end
    transitions = transitions .- minimum(transitions)
    chi_bar(T) = sum([exp(-transitions[i]/T) for i in 1:2:length(transitions)-1])/sum([exp(-transitions[i]/T) for i in 2:2:length(transitions)])
    χ_0 = target_acceptance_rate
    T_0 = 1.0
    try
        while abs(chi_bar(T_0) - χ_0) > 0.000001
            println(T_0)
            T_0 = T_0 * (log(chi_bar(T_0)) / log(χ_0 ))
        end
    catch 
        println("No energy decreasing transitions found!")
    end

    if isnan(T_0)
        println("NaN initial energy!")
        return NaN
    else
        return T_0
    end
end

