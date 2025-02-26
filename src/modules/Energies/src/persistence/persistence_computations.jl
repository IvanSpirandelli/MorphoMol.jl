function get_total_persistence_summed(dim_dgms::Vector{Matrix{Float64}}, weights::Vector{Float64} = [1.0, -1.0, 1.0])
    sum([get_total_persistence(dgm, weight) for (dgm, weight) in zip(dim_dgms, weights)])
end

function get_total_persistence(dgm, weight::Float64 = 1.0)
    if length(dgm) == 0
        return 0.0
    end
    weight * sum((dgm[:,2] - dgm[:,1]))
end

function get_death_by_birth_persistence(dgm, weight::Float64 = 1.0)
    if length(dgm) == 0
        return 0.0
    end
    weight * sum([dgm[i,2] / dgm[i,1] for i in 1:size(dgm)[1]])
end

function get_death_by_birth_persistence_summed(dim_dgms::Vector{Matrix{Float64}}, weights::Vector{Float64} = [1.0, -1.0, 1.0])
    sum([get_death_by_birth_persistence(dgm, weight) for (dgm, weight) in zip(dim_dgms, weights)])
end