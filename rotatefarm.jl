"""
    rotatefarm(coordinates, θ)

Rotates a wind farm based on the wind direction.

INPUTS
coordinates: 3 by n array of the n turbine coordinates
θ: the direction of wind source, measured in degrees clockwise from north

OUTPUTS
modified coordinates

"""


function rotatefarm(coordinates, θ)

    return [cos(-pi/180*(90+θ)) -sin(-pi/180*(90+θ)) 0; sin(-pi/180*(90+θ)) cos(-pi/180*(90+θ)) 0; 0 0 1]*coordinates

end
