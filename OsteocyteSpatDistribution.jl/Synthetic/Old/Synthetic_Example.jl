import OsteocyteSpatDistribution as OSD
using CairoMakie
using Random

Random.seed!(44)  # For reproducibility
# Example usage:
N = 10  # number of cells
M = 5    # number of radii
radii = range(10.0, 50.0, length=M)
a = 3.0 .+ 2.0 * rand(N)
b = 1.0 .+ 1.0 * rand(N)
cells = OSD.generate_ellipses(N, M, collect(radii), a, b)
ordered_cells = OSD.sort_cells(cells)

intersecting_rays = OSD.store_intersecting_rays(ordered_cells, 50, minimum(radii))
ray_lengths = OSD.calc_ray_lengths(intersecting_rays)
ray_angles = OSD.calc_ray_angles(ordered_cells, intersecting_rays)
forward_rays, backward_rays = OSD.separate_forward_backward_rays(ordered_cells, intersecting_rays)

# plotting code
fig = Figure(size = (800, 400), fontsize=22)

# Panel 1: Spatial distribution
ax = Axis(
    fig[1, 1],
    aspect = 1,
    title = "Synthetic Osteocytes",
    xlabel = "x (μm)",
    ylabel = "y (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)

# Plot concentric circles (radii)
for r in radii
    θ = range(0, 2π, length=200)
    x = r * cos.(θ)
    y = r * sin.(θ)
    lines!(ax, x, y, color=:gray, alpha=0.3, linewidth=1.5)
end

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
end

# Plot intersecting rays
for rays_from_cell in forward_rays
    for ray in rays_from_cell
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax, x, y, color=:crimson, alpha=0.7, linewidth=1.5)
    end
end

for rays_from_cell in backward_rays
    for ray in rays_from_cell
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax, x, y, color=:green, alpha=0.7, linewidth=1.5)
    end
end

# Panel 2: Ray length histogram
ax2 = Axis(
    fig[1, 2],
    title = "Distribution of Ray Lengths",
    xlabel = "Ray Length (μm)",
    ylabel = "Count",
    xticklabelsize = 18,
    yticklabelsize = 18,
    backgroundcolor = :white
)
hist!(ax2, ray_lengths, bins=20, color=:crimson, strokewidth=1.5, strokecolor=:black)
# Remove top/right spines for a cleaner look
hidespines!(ax, :t, :r)
hidespines!(ax2, :t, :r)

fig

"""

fig2 = Figure(size = (600, 600), fontsize=22)

intersecting_rays_10 = OSD.store_intersecting_rays(ordered_cells, 10, minimum(radii))
intersecting_rays_20 = OSD.store_intersecting_rays(ordered_cells, 20, minimum(radii))
intersecting_rays_30 = OSD.store_intersecting_rays(ordered_cells, 30, minimum(radii))
intersecting_rays_40 = OSD.store_intersecting_rays(ordered_cells, 40, minimum(radii))
intersecting_rays_50 = OSD.store_intersecting_rays(ordered_cells, 50, minimum(radii))

# Panel 3: Rays per cells
ax = Axis(
    fig2[1, 1],
    title = "Rays per Cell",
    xlabel = "Osteocyte",
    ylabel = "Count",
    xticklabelsize = 18,
    yticklabelsize = 18,
    backgroundcolor = :white
)
rays_per_cell_10 = [length(rays) for rays in intersecting_rays_10]
rays_per_cell_20 = [length(rays) for rays in intersecting_rays_20]
rays_per_cell_30 = [length(rays) for rays in intersecting_rays_30]
rays_per_cell_40 = [length(rays) for rays in intersecting_rays_40]
rays_per_cell_50 = [length(rays) for rays in intersecting_rays_50]

lines!(ax, 1:N, rays_per_cell_10, color=:dodgerblue, label="10 μm Rays", linewidth = 2)
lines!(ax, 1:N, rays_per_cell_20, color=:crimson, label="20 μm Rays", linewidth = 2)
lines!(ax, 1:N, rays_per_cell_30, color=:black, label="30 μm Rays", linewidth = 2)
lines!(ax, 1:N, rays_per_cell_40, color=:green, label="40 μm Rays", linewidth = 2)
lines!(ax, 1:N, rays_per_cell_50, color=:purple, label="50 μm Rays", linewidth = 2)

axislegend(position = :lt)
hidespines!(ax, :t, :r)

fig2


# Animation: highlight each cell in sequence
frames = 60
cell_indices = 1:N

record(fig, "cell_animation.mp4", cell_indices; framerate=1) do i
    scatter!(ax, [ordered_cells[i].center[1]], [ordered_cells[i].center[2]], color=:orange, markersize=30)
    sleep(0.05)  # Optional: slow down for visibility
    if i > 1
        scatter!(ax, [ordered_cells[i-1].center[1]], [ordered_cells[i-1].center[2]], color=:white, markersize=30)
    end
end
"""

# plotting code
fig2 = Figure(size = (800, 400), fontsize=22)

# Panel 1: Spatial distribution
ax = Axis(
    fig2[1, 1],
    aspect = 1,
    title = "Forward Rays",
    xlabel = "x (μm)",
    ylabel = "y (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)

# Plot concentric circles (radii)
for r in radii
    θ = range(0, 2π, length=200)
    x = r * cos.(θ)
    y = r * sin.(θ)
    lines!(ax, x, y, color=:gray, alpha=0.3, linewidth=1.5)
end

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
    lines!(ax, x_rot, y_rot, alpha=0.8, linewidth=2)
end

# Plot intersecting rays
for rays_from_cell in forward_rays
    for ray in rays_from_cell
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax, x, y, color=:crimson, alpha=0.7, linewidth=1.5)
    end
end

ax2 = Axis(
    fig2[1, 2],
    aspect = 1,
    title = "Backward Rays",
    xlabel = "x (μm)",
    ylabel = "y (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)

# Plot concentric circles (radii)
for r in radii
    θ = range(0, 2π, length=200)
    x = r * cos.(θ)
    y = r * sin.(θ)
    lines!(ax2, x, y, color=:gray, alpha=0.3, linewidth=1.5)
end

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
    lines!(ax2, x_rot, y_rot, alpha=0.8, linewidth=2)
end

for rays_from_cell in backward_rays
    for ray in rays_from_cell
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax2, x, y, color=:green, alpha=0.7, linewidth=1.5)
    end
end

fig2

intersection_counts = OSD.count_ray_intersections(cells, backward_rays)