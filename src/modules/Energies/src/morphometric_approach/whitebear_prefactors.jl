function get_prefactors(rs::Float64, η::Float64)
    [ 
    pressure(rs, η), # Volume
    sigma(rs, η), # Area
    kappa(rs, η), # Mean curvature 
    kappa_bar(η), # Gaussian curvature
    ]
end

# These prefactors are equivalent to doubling the negative contributions to the mean curvature 
# for a union of hard spheres of radius r. (This was a bug in a previous version of AlphaMol) 
function get_bugged_prefactors(r::Float64, rs::Float64, η::Float64)
    pf_wb = get_prefactors(rs, η)
    [pf_wb[1], pf_wb[2] - (pf_wb[3]/(r + rs)), 2.0*pf_wb[3], pf_wb[4]]
end

""" Prefactors.
    Given in equation.24 of Density functional theory for hard-sphere mixtures:
    the White Bear version mark II.
"""

function pressure(rs::Float64, η::Float64)
    fraction = (1+η+η^2-η^3)/(1-η)^3
    scaling = η*3/4/pi/rs^3
    scaling * fraction
end

function sigma(rs::Float64, η::Float64)
    fraction_1 = (1+2*η + 8*η^2 - 5*η^3)/(3*(1-η)^3)
    fraction_2 = (log(1-η)/(3*η))
    scaling = -η*3/4/pi/rs^2 #Note the MINUS
    scaling * (fraction_1 + fraction_2)
end

function kappa(rs::Float64, η::Float64)
    fraction_1 = (4-10*η+20*η^2-8*η^3)/(3*(1-η)^3)
    fraction_2 = (4*log(1-η)/(3*η))
    scaling = η*3/4/pi/rs
    scaling * (fraction_1 + fraction_2)
end

function kappa_bar(η::Float64)
    fraction_1 = (-4+11*η-13*η^2+4*η^3)/(3*(1-η)^3)
    fraction_2 = (4*log(1-η))/(3*η)
    scaling = η*3/4/pi
    scaling * (fraction_1 - fraction_2)
end
