function state_to_poly(flat_realization::Vector, radii::Vector, filepath::String, n_mol::Int, n_atoms_per_mol::Int)
    open(string(filepath, ".poly"), "w") do io
        println(io,"POINTS")
        color = ""
        for i in 0:n_mol-1
            color = "c($(rand()),$(rand()),$(rand())"
            for j in 0:n_atoms_per_mol-1
                atom_id = (i * n_atoms_per_mol + j)
                x = atom_id*3 + 1
                y = atom_id*3 + 2
                z = atom_id*3 + 3
                println(io, "$(i * n_atoms_per_mol + j + 1): $(flat_realization[x]) $(flat_realization[y]) $(flat_realization[z]) $(radii[atom_id+1]) $(color), $(radii[atom_id+1]))")
            end
        end
        println(io,"POLYS")
        println(io,"END")
    end
end

function poly_to_state(filepath)
    coordinate_data = readlines(filepath)[2:end-2]
    flat_realization = Vector{Float64}([])
    radii = Vector{Float64}([])
    # Assuming lines are formatted like "i: x y z r whatever"
    for line in coordinate_data
        line_array = split(line, " ")
        push!(flat_realization, parse(Float64, line_array[2]))
        push!(flat_realization, parse(Float64, line_array[3]))
        push!(flat_realization, parse(Float64, line_array[4]))
        push!(radii, parse(Float64, line_array[5]))
    end
    flat_realization, radii
end
