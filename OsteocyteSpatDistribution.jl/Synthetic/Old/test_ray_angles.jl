import OsteocyteSpatDistribution as OSD
using CairoMakie
using Random

cell1 = OSD.Ellipse([0.0, 0.0], [2.0, 0.5], 0.0)
cell2 = OSD.Ellipse([0.0, 5.0], [2.0, 0.5], 0.0)
cell3 = OSD.Ellipse([20.0, 2.5], [2.0, 0.5], 0.0)

cells = [cell1, cell2, cell3]
all_rays = OSD.generate_rays_from_cell(cell1, 10, 100)
intersecting_rays = OSD.store_intersecting_rays(cells, 10, 0.0)

ray_lengths = OSD.calc_ray_lengths(intersecting_rays)
ray_angles = OSD.calc_ray_angles(cells, intersecting_rays)

fig = Figure(size = (500, 500), fontsize=22)
ax = Axis(
    fig[1, 1],
    title = "Synthetic Osteocytes",
    xlabel = "x (μm)",
    ylabel = "y (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)

# Plot ellipses (cells)
for e in cells
    θ = range(0, 2π, length=200)
    a, b = e.axes
    xs = a * cos.(θ)
    ys = b * sin.(θ)
    R = [cos(e.angle) -sin(e.angle); sin(e.angle) cos(e.angle)]
    pts = [R * [xs[i], ys[i]] + e.center for i in 1:length(θ)]
    x_rot = [p[1] for p in pts]
    y_rot = [p[2] for p in pts]
    lines!(ax, x_rot, y_rot, color=:dodgerblue, alpha=0.8, linewidth=2)
    #scatter!(ax, x_rot[10], y_rot[10], color=:green, markersize = 15)
end


for ray in all_rays
    x = [ray[1][1], ray[2][1]]
    y = [ray[1][2], ray[2][2]]
    lines!(ax, x, y, color=:grey, alpha=0.5, linewidth=1.5)
end


# Plot intersecting rays
for rays_from_cell in intersecting_rays
    for ray in rays_from_cell
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax, x, y, color=:red, alpha=0.5, linewidth=1.5)
    end
end

fig

intersection_counts = OSD.count_ray_intersections(cells, intersecting_rays)
