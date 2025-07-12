import OsteocyteSpatDistribution as OSD
using CairoMakie
using Random

Random.seed!(45)  # For reproducibility
# Example usage:
N = 20  # number of cells
M = 5    # number of radii
radii = range(10.0, 50.0, length=M)
a = 3.0 .+ 2.0 * rand(20)
b = 1.0 .+ 1.0 * rand(20)
cells = OSD.generate_ellipses(N, M, collect(radii), a, b)

intersecting_rays = OSD.store_intersecting_rays(cells, 50, minimum(radii))
ray_lengths = OSD.calc_ray_lengths(intersecting_rays)


# plotting code

fig = Figure(size = (1000, 500), fontsize=22)

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
for rays_from_cell in intersecting_rays
    for ray in rays_from_cell
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax, x, y, color=:crimson, alpha=0.7, linewidth=1.5)
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
