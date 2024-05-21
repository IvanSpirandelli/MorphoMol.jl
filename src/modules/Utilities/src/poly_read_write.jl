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
