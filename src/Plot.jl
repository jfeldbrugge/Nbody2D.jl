" Mesh plot of the N-body simulation"
function plotMesh!(state, box, color = :black, margin = 0., title = "")
	x = state.position[:,:,1]
	y = state.position[:,:,2]
	
	Plots.plot!(x, y, c = color, aspect_ratio=:equal, 
    xlims=(-margin * box.L, (1. + margin) * box.L),
    ylims=(-margin * box.L, (1. + margin) * box.L), 
    label=false, size=(500, 500), linewidth=0.1, title = title)
	
	Plots.plot!(x', y', c = color, aspect_ratio=:equal, 
    xlims=(-margin * box.L, (1. + margin) * box.L),
    ylims=(-margin * box.L, (1. + margin) * box.L), 
    label=false, linewidth=0.1)
end

" Mesh plot of the N-body simulation"
function plotMesh(state, box, color = :black, margin = 0., title = "")
    plot()
    return plotMesh!(state, box, color, margin, title)
end
" Plot Caustic Skeleton"
function plotSkeleton!(skeleton::CausticSkeleton, plotA2 = true)
    if plotA2
        plot!([(sim[:,1], sim[:,2]) for sim in skeleton.A2], label = false, color = :blue, linewidth = 2)
    end

    plot!([(sim[:,1], sim[:,2]) for sim in skeleton.A3], label = false, color = :red, linewidth = 2)
    plot!()    
end

" Plot Caustic Skeleton"
function plotSkeleton(skeleton::CausticSkeleton, plotA2 = true)
    plot(aspect_ratio=:equal)
    plotSkeleton!(skeleton, plotA2)
end
