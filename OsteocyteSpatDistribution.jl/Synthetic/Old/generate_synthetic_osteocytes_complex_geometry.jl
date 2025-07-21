using Interpolations

# Suppose curve_points is a Nx2 array of (x, y) points along the closed curve
curve_points = [cos(θ), sin(θ) for θ in range(0, 2π, length=100)] # Example: unit circle with 100 points
curve_points = hcat(curve_points...)' # Convert to a Nx2 matrix
t = range(0, 1; length=size(curve_points, 1)+1)[1:end-1] # parameterization

# Interpolate the curve
itp_x = LinearInterpolation(t, curve_points[:,1], extrapolation_bc=Periodic())
itp_y = LinearInterpolation(t, curve_points[:,2], extrapolation_bc=Periodic())

function curve_func(s)
    [itp_x(s), itp_y(s)]
end

function curve_tangent(s, δ=1e-5)
    (curve_func(s+δ) - curve_func(s-δ)) / (2δ)
end

function generate_ellipses_on_curve(N, curve_func, curve_tangent; a=2.0, b=1.0, max_attempts=100)
    ellipses = Ellipse[]
    placed = 0
    attempts = 0
    while placed < N && attempts < max_attempts
        s = rand() # random parameter along curve
        center = curve_func(s)
        tangent = curve_tangent(s)
        angle = atan(tangent[2], tangent[1]) # orientation along tangent
        new_ellipse = Ellipse(center, [a, b], angle)
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



# Helper to evaluate the curve at parameter t in [0,1]
function eval_curve(curve, t)
    println(curve)
    println(typeof(curve))
    println(all(x -> isa(x, Tuple{Float64,Float64}), curve))
    if typeof(curve) <: Function
        return curve(t)
    elseif all(x -> isa(x, Tuple{Float64,Float64}), curve)
        n = length(curve)
        tt = clamp(t * (n-1) + 1, 1, n)
        i = floor(Int, tt)
        j = clamp(i+1, 1, n)
        α = tt - i
        x = (1-α)*curve[i][1] + α*curve[j][1]
        y = (1-α)*curve[i][2] + α*curve[j][2]
        return (x, y)
    end
end