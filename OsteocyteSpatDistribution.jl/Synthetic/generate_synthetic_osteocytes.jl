using Random
using LinearAlgebra
using CairoMakie

struct Ellipse
    center::Vector{Float64}
    axes::Vector{Float64}  # [a, b] semi-axes lengths
    angle::Float64         # orientation angle in radians
end

"""
    generate_ellipses(N, M, radii; a=2.0, b=1.0)

Generate N ellipses randomly distributed at M radii (given as a vector `radii`).
Each ellipse is placed in the direction of the radius, with random orientation.
Returns a vector of Ellipse structs.
"""
function generate_ellipses(N, M, radii; a=2.0, b=1.0, max_attempts=100)
    ellipses = Ellipse[]
    placed = 0
    attempts = 0
    while placed < N && attempts < max_attempts
        # Choose a radius
        r = radii[rand(1:M)]
        # Random direction on the circle
        θ = 2π * rand()
        x = r * cos(θ)
        y = r * sin(θ)
        center = [x, y]
        # Orientation angle is the direction from origin to center
        angle = θ + pi/2
        new_ellipse = Ellipse(center, [a, b], angle)
        # Check for overlap (using bounding circles for simplicity)
        min_dist = 2 * maximum([a, b])
        overlap = any(norm(e.center - center) < min_dist for e in ellipses)
        if !overlap
            push!(ellipses, new_ellipse)
            placed += 1
        end
        attempts += 1
    end
    if placed < N
        @warn "Could only place $placed out of $N ellipses without overlap."
    end
    return ellipses
end

# Example usage:
N = 20  # number of cells
M = 5    # number of radii
radii = range(10.0, 50.0, length=M)
cells = generate_ellipses(N, M, collect(radii))


fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], aspect = 1)

# Plot circles at each radius
for r in radii
    θ = range(0, 2π, length=100)
    x = r * cos.(θ)
    y = r * sin.(θ)
    lines!(ax, x, y, color=:gray, alpha=0.2)
end

# Plot ellipses (cells)
for e in cells
    θ = range(0, 2π, length=100)
    a, b = e.axes
    # Parametric ellipse
    xs = a * cos.(θ)
    ys = b * sin.(θ)
    # Rotate and translate
    R = [cos(e.angle) -sin(e.angle); sin(e.angle) cos(e.angle)]
    pts = [R * [xs[i], ys[i]] + e.center for i in 1:length(θ)]
    x_rot = [p[1] for p in pts]
    y_rot = [p[2] for p in pts]
    lines!(ax, x_rot, y_rot, color=:blue, alpha=0.5)
end

"""
for cell in cells
    rays = generate_rays_from_cell(cell, 50)
    for ray in rays
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax, x, y, color=:grey, alpha=0.5)
    end
end
"""

intersecting_rays = store_intersecting_rays(cells, 50, minimum(radii))

for rays_from_cell in intersecting_rays
    for ray in rays_from_cell
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax, x, y, color=:red, alpha=0.5)
    end
end

fig

