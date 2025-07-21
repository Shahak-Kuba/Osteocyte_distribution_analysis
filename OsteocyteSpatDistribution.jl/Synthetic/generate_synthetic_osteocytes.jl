
function generate_osteocyte_points(N::Int, R::AbstractVector{<:Real}, distribution::String)
    curves = [t -> (r * cos(2π * t), r * sin(2π * t)) for r in R]
    points = Vector{Tuple{Float64, Float64}}()
    distances = Float64[]
    tangents = Vector{Tuple{Float64, Float64}}()
    normals = Vector{Tuple{Float64, Float64}}()
    if distribution == "Uniform"
        n_curves = length(curves)
        n_per_curve = fill(div(N, n_curves), n_curves)
        for i in 1:rem(N, n_curves)
            n_per_curve[i] += 1
        end
        for (curve_idx, (curve, n_pts)) in enumerate(zip(curves, n_per_curve))
            r = R[curve_idx]
            for i in 1:n_pts
                t = (i-1) / n_pts
                x, y = curve(t)
                push!(points, (x, y))
                push!(distances, r)
                # Tangent: derivative of (r*cos(2πt), r*sin(2πt)) w.r.t t
                tx = -r * 2π * sin(2π * t)
                ty =  r * 2π * cos(2π * t)
                norm_tan = sqrt(tx^2 + ty^2)
                push!(tangents, (tx / norm_tan, ty / norm_tan))
                # Inward normal: for a circle, it's just (-x/r, -y/r)
                push!(normals, (-x / r, -y / r))
            end
        end
    elseif distribution == "Random"
        for _ in 1:N
            curve_idx = rand(1:length(curves))
            curve = curves[curve_idx]
            r = R[curve_idx]
            t = rand()
            x, y = curve(t)
            push!(points, (x, y))
            push!(distances, r)
            tx = -r * 2π * sin(2π * t)
            ty =  r * 2π * cos(2π * t)
            norm_tan = sqrt(tx^2 + ty^2)
            push!(tangents, (tx / norm_tan, ty / norm_tan))
            push!(normals, (-x / r, -y / r))
        end
    else
        throw(ArgumentError("Unknown distribution: $distribution"))
    end
    # Create Osteocyte struct for each point
    Osteocytes = Vector{Osteocyte}()
    for (point, dist, tan, norm) in zip(points, distances, tangents, normals)
        # Create Osteocyte struct for each point
        osteocyte = Osteocyte(point, dist, tan, norm)
        push!(Osteocytes, osteocyte)
        # Here you can store or process the osteocyte as needed
    end
    return Osteocytes
end

N = 30; R = range(10.0, 50.0, length=5); distribution = "Random"
points, distances, tangents, normals = generate_osteocyte_points(N,R, distribution)


using CairoMakie

f = Figure(size = (400, 400), fontsize=22)

# --- Left plot: Circle points ---
ax = Axis(
    f[1, 1],
    aspect = 1,
    title = "Circle (Uniform)",
    xlabel = "x (μm)",
    ylabel = "y (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
# Plot concentric circles (radii)
for r in R
    θ = range(0, 2π, length=200)
    x = r * cos.(θ)
    y = r * sin.(θ)
    lines!(ax, x, y, color=:gray, alpha=0.3, linewidth=1.5)
end 
# Plot random points on the circle
scatter!(ax, points, color=:dodgerblue, markersize=10)
f

