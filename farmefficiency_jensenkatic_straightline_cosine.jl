"""

    farmefficiency_jensenkatic_straightline_cosine(coordinates, D, a, α, k)

Calculate the efficiency of a given farm with a given wind direction, using straight-line approximations for turbines partially in wakes, and using Jensen's cosine method for the velocity variation across a wake.

# Arguments
- 'coordinates::Array{Float64,2}': 2 by n array with the coordinates for each wind turbine
- 'D::Float64': diameter of a turbine
- 'a::Float64': initial velocity deficit
- 'α::Float64': entrainment coefficient
- 'k::Float64': amount of lateral spread of wake per unit longitudinal distance

"""

function farmefficiency_jensenkatic_straightline_cosine(coordinates, D, a, α, k)


    function areafractioninwake_cosine(x, wake, turbine)

        function f(y, w, D, y_offset)
            if abs(y) <= w
                return ((1 + cos(pi*y/w))/2)^(.8) * 2*sqrt(maximum(((D/2)^2 - (y - y_offset)^2, 0)))
            else
                return 0
            end
        end

        # bounds are max of y-D/2 or -D/2
        if (turbine[1] < wake[1] && turbine[2] < wake[1]) || (turbine[1] > wake[2] && turbine[2] > wake[2])
            return 0
        else
            lb = maximum((wake[1], turbine[1]))
            ub = minimum((wake[2], turbine[2]))
            numsegments = 1000  # number of segments
            g, y = zeros(numsegments+1), zeros(numsegments+1)
            y[1] = lb
            g[1] = f(y[1], (wake[2]-wake[1])/2, turbine[2]-turbine[1], (turbine[2]+turbine[1])/2-(wake[2]+wake[1])/2)
            velocitydef_fraction = zeros(numsegments)
            for i = 1:numsegments
                y[i+1] = lb + i*(ub-lb)/numsegments
                g[i+1] = f(y[i+1], (wake[2]-wake[1])/2, turbine[2]-turbine[1], (turbine[2]+turbine[1])/2-(wake[2]+wake[1])/2)
                velocitydef_fraction[i] = (g[i] + g[i+1])/2 * (ub-lb)/numsegments / (pi/4 * (turbine[2]-turbine[1])^2)
            end
            return sum(velocitydef_fraction) 
        end
    end

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
                    if !((turbine_lowerbound < wake_lowerbound && turbine_upperbound < wake_lowerbound) || (turbine_lowerbound > wake_upperbound && turbine_upperbound > wake_upperbound))
                        areafractioninwake = areafractioninwake_cosine(connectvec[1], [wake_lowerbound; wake_upperbound], [turbine_lowerbound; turbine_upperbound])
                        energydeftotal += areafractioninwake*singlewakevelocitydeficit(connectvec[1], D, a, α)^2
                    end

                end

            end

        end

        wakevelocity[i] = 1 - sqrt(energydeftotal)
        energydeftotal = 0

    end


    return sum(wakevelocity.^3) / n

end
