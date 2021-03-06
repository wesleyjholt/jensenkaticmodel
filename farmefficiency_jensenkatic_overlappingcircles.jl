"""

    farmefficiency_jensenkatic_overlappingcircles(coordinates, D, a, α, k)

Calculate the efficiency of a given farm with a given wind direction, using overlapping circles for turbines partially in wakes.

# Arguments
- 'coordinates::Array{Float64,2}': 2 by n array with the coordinates for each wind turbine
- 'D::Float64': diameter of a turbine
- 'a::Float64': initial velocity deficit
- 'α::Float64': entrainment coefficient
- 'k::Float64': amount of lateral spread of wake per unit longitudinal distance
    
"""

function farmefficiency_jensenkatic_overlappingcircles(coordinates, D, a, α, k)

include("singlewakevelocitydeficit.jl")
include("intersecting_circles_area.jl")

    # energy deficit
    energydeftotal = 0
    n = length(coordinates[1,:])
    wakevelocity = zeros(n,1)

    for i = 1:n
        for j = 1:n

            if i !== j  # effect of a turbine on itself is zero
                connectvec = coordinates[:,i] - coordinates[:,j] # position vector from the center of turbine j to i
                wake_upperbound = D/2 + k*connectvec[1]          # y-position of wake upper bound w.r.t. center of turbine j (at x-position of turbine i)
                wake_lowerbound = -D/2 - k*connectvec[1]         # y-position of wake lower bound w.r.t. center of turbine j (at x-position of turbine i)
                turbine_upperbound = connectvec[2] + D/2         # y-position of turbine i upper bound w.r.t. center of turbine j
                turbine_lowerbound = connectvec[2] - D/2         # y-position of turbine i lower bound w.r.t. center of turbine j

                if connectvec[1] > 0
                     # turbine i is downstream from turbine j

                    if wake_lowerbound < turbine_upperbound < wake_upperbound && wake_lowerbound < turbine_lowerbound < wake_upperbound
                        # turbine i is completely in the wake of turbine j
                        energydeftotal += singlewakevelocitydeficit(connectvec[1], D, a, α)^2

                    elseif wake_lowerbound < turbine_upperbound < wake_upperbound && turbine_lowerbound < wake_lowerbound
                        # turbine i is partially in the wake of turbine j
                        # (lower portion of turbine i is outside the wake)

                            areafractioninwake = intersecting_circles_area(D/2, wake_upperbound, abs(connectvec[2])) / (pi * (D/2)^2)
                            energydeftotal += areafractioninwake*singlewakevelocitydeficit(connectvec[1], D, a, α)^2

                    elseif turbine_upperbound > wake_upperbound && wake_lowerbound < turbine_lowerbound < wake_upperbound
                        # turbine i is partially in the wake of turbine j
                        # (upper portion of turbine i is outside the wake)

                            areafractioninwake = intersecting_circles_area(D/2, wake_upperbound, abs(connectvec[2])) / (pi * (D/2)^2)
                            energydeftotal += areafractioninwake*singlewakevelocitydeficit(connectvec[1], D, a, α)^2
                            # if areafractioninwake > 1
                            #     error("Area fraction in wake is greater than one")
                            # end

                    end

                end

            end

        end

        wakevelocity[i] = 1 - sqrt(energydeftotal)
        if wakevelocity[i] < 0
            error("Wake velocity deficit is negative")
        end
        energydeftotal = 0

    end


    return sum(wakevelocity.^3) / n

end
