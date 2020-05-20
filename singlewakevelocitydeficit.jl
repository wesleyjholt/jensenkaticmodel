"""

        singlewakevelocitydeficit(x, D, a, alphavar)

Calculates the wind velocity at a turbine relative to the ambient wind velocity.

INPUTS
x:      distance between turbines (vector, each entry for the effect from a different turbine)
D:      rotor diameter
a:      initial velocity deficit
alphavar:  entrainment coefficient

"""

# singlewakevelocitydeficit
# Author: Wesley Holt
# Date: 5/14/20

function singlewakevelocitydeficit(x, y, D, a, α, k)

        θ_offaxis = abs(atan(y,x))
        θ_max = atan(k)
        return 2*a./(1 .+ 2*α.*x./D).^2 #* (1 + cos(180/θ_max * θ_offaxis * pi/180))/2

end
