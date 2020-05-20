"""
    intersecting_circles_area(r1, r2, d)

Compute the area of two overlapping circles.

# Arguments
- 'r1::Float': the radius of the first circle
- 'r2::Float': the radius of the second circle
- 'd::Float': the distance between the centers of the two circles

"""

function intersecting_circles_area(r1, r2, d)

    if d < r1 + r2
        return r1^2 * acos((d^2 + r1^2 - r2^2)/(2*d*r1)) + r2^2 * acos((d^2 + r2^2 - r1^2)/(2*d*r2)) - 1/2 * sqrt((-d + r1 + r2)*(d + r1 - r2)*(d - r1 + r2)*(d + r1 + r2))
    else
        return 0
    end
    
end