# katicfig4&5
# Author: Wesley Holt
# Date: 5/20/20

using Plots
using DataFrames
include("singlewakevelocitydeficit.jl")

# INITIAL PARAMETERS
    D = 40      # diameter of one turbine
    u = 9.3       # ambient wind speed
    a = 1/3     
    α1 = 0.075     
    α2 = 0.11      
    k = 0.1     #0.06

# VELOCITY AT THE MASTS
    v1 = u * (1 - singlewakevelocitydeficit(2.5*D, D, a, α1) )
    v2 = u * (1 - singlewakevelocitydeficit(4*D, D, a, α1) )
    v3a = u * (1 - sqrt( singlewakevelocitydeficit(6*D, D, a, α1)^2 + singlewakevelocitydeficit(D, D, a, α1)^2 ) )
    v3b = u * (1 - sqrt( singlewakevelocitydeficit(6*D, D, a, α1)^2 + singlewakevelocitydeficit(D, D, a, α2)^2 ) )
    v4a = u * (1 - sqrt( singlewakevelocitydeficit(7.5*D, D, a, α1)^2 + singlewakevelocitydeficit(2.5*D, D, a, α1)^2 ) )
    v4b = u * (1 - sqrt( singlewakevelocitydeficit(7.5*D, D, a, α1)^2 + singlewakevelocitydeficit(2.5*D, D, a, α2)^2 ) )

# COMPARE MY RESULTS WITH KATIC'S RESULTS
    three_a = 3.255595468361427
    three_b = 4.592889380123423
    four_a = 5.38786036658377
    one_a = 6.321820023947684
    four_b = 6.321820023947684
    two_a = 7.25550336188634

    
    scatter([one_a; two_a; three_a; three_b; four_a; four_b], 45*ones(6),
            xlim = (2,8),
            ylim = (30,60))
    scatter!([v1, v2, v3a, v3b, v4a, v4b], 45*ones(6))

    results = DataFrame(Mast = ["One","Two","Three","Three","Four","Four"],
                        α = [0.075,0.075,0.075,0.11,0.075,0.11,],
                        Katic_wind_velocity = [one_a, two_a, three_a, three_b, four_a, four_b],
                        My_wind_velocity = [v1, v2, v3a, v3b, v4a, v4b],
                        Percent_error = [100*(one_a-v1)/one_a, 100*(two_a-v2)/two_a, 100*(three_a-v3a)/three_a, 100*(three_b-v3b)/three_b, 100*(four_a-v4a)/four_a, 100*(four_b-v4b)/four_b]
                        )
    display(results)