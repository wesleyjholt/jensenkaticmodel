"""
    rotatefarm(coordinates, θ)

Rotates a wind farm based on the wind direction.

# Arguments
- 'coordinates::Array{Float64,2}': 3 by n array of the n turbine coordinates
- 'θ::Float64': the direction of wind source, measured in degrees clockwise from north

"""

function rotatefarm(coordinates, θ)

    return [cos(-pi/180*(90+θ)) -sin(-pi/180*(90+θ)) 0; sin(-pi/180*(90+θ)) cos(-pi/180*(90+θ)) 0; 0 0 1]*coordinates

end
