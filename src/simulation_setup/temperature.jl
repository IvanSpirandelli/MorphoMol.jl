# This function takes a sequence of energy evaluations and a a target acceptance rate 
# to compute the temperature needed to achieve the desired targe rate in another simulation.
# This requires other paramters in the simulation to be the same as the original simulation.
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