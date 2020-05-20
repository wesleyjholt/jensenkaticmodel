# katicfig3_jensenkaticmodel
# Author: Wesley Holt
# Date: 5/14/20

using Plots, LinearAlgebra, Statistics, CSV, DataFrames
include("farmefficiency_jensenkatic_overlappingcircles.jl")
include("farmefficiency_jensenkatic_straightline.jl")
include("farmefficiency_jensenkatic_straightline_cosine.jl")
include("rotatefarm.jl")

#---INITIAL PARAMETERS---

    D = 20      # diameter of one turbine
    h = D       # turbine height
    u = 1       # ambient wind speed
    a = 1/3     
    α = 0.0775    #0.08   
    k = 0.05    #0.06

    coordinates = convert(Matrix, CSV.read("katicfig3_coordinates.csv", header=false))

#---FIND THE EFFICIENCY OVER MULTIPLE WIND DIRECTIONS---

    θ_sweeplength = 200
    θ_sweep = range(0, stop=45, length = θ_sweeplength)
    efficiency_overlappingcircles, efficiency_straightline, efficiency_straightline_cosine = zeros(θ_sweeplength), zeros(θ_sweeplength), zeros(θ_sweeplength)

    for i = 1:θ_sweeplength
        θ = θ_sweep[i]
        # call farmefficiency function
        efficiency_overlappingcircles[i] = farmefficiency_jensenkatic_overlappingcircles(rotatefarm(coordinates, θ), D, a, α, k)
        efficiency_straightline[i] = farmefficiency_jensenkatic_straightline(rotatefarm(coordinates, θ), D, a, α, k)
        efficiency_straightline_cosine[i] = farmefficiency_jensenkatic_straightline_cosine(rotatefarm(coordinates, θ), D, a, α, k)
    end

#--- DISPLAY RESULTS---

        font_size = 10
        line_width = 3

        # overlapping circles plot
        katicfig3_plot = plot(θ_sweep, efficiency_overlappingcircles,
            xlabel = "Wind direction (degrees)",
            ylabel = "E/Eo",
            linewidth = line_width,
            legend = false,
            label = "My Results",
            xtickfont=font(font_size -2),
            ytickfont=font(font_size -2),
            guidefont=font(font_size + 2),
            legendfont=font(font_size),
            xlim = [0, 45],
            ylim = [.5, 1],
            xticks = [0,20,40],
            yticks = [.5,1],
            grid = false)

        # straight line plot
        katicfig3_plot = plot!(θ_sweep, efficiency_straightline,
            linewidth = line_width,
            label = "My Results, straight line"  )

        # straight line with cosine model
        katicfig3_plot = plot!(θ_sweep, efficiency_straightline_cosine,
            linewidth = line_width,
            label = "My Results, straight line with cosine"  )

        # Katic's results    
        katicfig3_data = CSV.read("katicfig3_data.csv", header=false)
        katicfig3_plot = plot!(katicfig3_data[:,1], katicfig3_data[:,2],
            linestyle = :dash,
            linewidth = line_width + 1,
            linecolor = :black,
            label = "Katic's Results"  )
        
        
        display(katicfig3_plot)

