using Random
using LinearAlgebra

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

function generate_ellipses(N, M, radii, a_vec, b_vec; max_attempts=100)
    @assert length(a_vec) == N "a_vec must have length N"
    @assert length(b_vec) == N "b_vec must have length N"
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
        a = a_vec[placed+1]
        b = b_vec[placed+1]
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

"""
    atan2pi(y, x)

Returns the angle in radians between the positive x-axis and the point (x, y),
in the range [0, 2π).
"""
function atan2pi(y, x)
    θ = atan(y, x)
    return θ < 0 ? θ + 2π : θ
end

"""
    sort_ellipses_by_radius_and_angle(ellipses)

Returns a new vector of ellipses sorted first by descending distance from the origin (furthest first),
then by the angle (in radians) each ellipse's center makes with the x-axis (increasing from 0 to 2π).
"""
function sort_cells(cells) 
    return sort(cells; by = e -> (-norm(e.center)))
end
