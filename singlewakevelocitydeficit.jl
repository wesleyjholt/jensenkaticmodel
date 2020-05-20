"""

        singlewakevelocitydeficit(x, D, a, α)

Calculates the wind velocity at a turbine relative to the ambient wind velocity.

# Arguments
- 'x::Float64': perpendicular-to-wind-direction distance between turbines
- 'D::Float64': rotor diameter
- 'a::Float64': initial velocity deficit
- 'α::Float64': entrainment coefficient

"""

# singlewakevelocitydeficit
# Author: Wesley Holt
# Date: 5/14/20

function singlewakevelocitydeficit(x, D, a, α)

        return 2*a./(1 .+ 2*α.*x./D).^2 #* (1 + cos(180/θ_max * θ_offaxis * pi/180))/2

end
