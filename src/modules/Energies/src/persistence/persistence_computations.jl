function get_total_persistence(dgm)
    if length(dgm) == 0
        return 0.0
    end
    sum((dgm[:,2] - dgm[:,1]))
end