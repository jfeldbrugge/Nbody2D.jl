" Mesh plot of the N-body simulation"
function plotMesh(state, box, color = :black, margin = 0.)
	x = state.position[:,:,1]
	y = state.position[:,:,2]
	
	plot(x, y, c = color, aspect_ratio=:equal, 
    xlims=(-margin * box.L, (1. + margin) * box.L),
    ylims=(-margin * box.L, (1. + margin) * box.L), 
    label=false, size=(500, 500), linewidth=0.1)
	
	plot!(x', y', c = color, aspect_ratio=:equal, 
    xlims=(-margin * box.L, (1. + margin) * box.L),
    ylims=(-margin * box.L, (1. + margin) * box.L), 
    label=false, linewidth=0.1)
end

" Mesh plot of the N-body simulation"
function plotMesh!(state, box, color = :black, margin = 0.)
	x = state.position[:,:,1]
	y = state.position[:,:,2]
	
	plot!(x, y, c = color, aspect_ratio=:equal, 
    xlims=(-margin * box.L, (1. + margin) * box.L),
    ylims=(-margin * box.L, (1. + margin) * box.L), 
    label=false, size=(500, 500), linewidth=0.1)
	
	plot!(x', y', c = color, aspect_ratio=:equal, 
    xlims=(-margin * box.L, (1. + margin) * box.L),
    ylims=(-margin * box.L, (1. + margin) * box.L), 
    label=false, linewidth=0.1)
end