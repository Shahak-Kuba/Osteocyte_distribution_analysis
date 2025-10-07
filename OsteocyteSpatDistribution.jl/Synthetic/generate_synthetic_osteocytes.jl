"""
    numerical_derivative(curve::Function, t::Real; h::Real=1e-5)

Numerically approximates the derivative of a parametric curve at parameter `t` using central differences.
Returns a tuple (dx/dt, dy/dt).
"""
function numerical_derivative(curve::Function, t::Real; h::Real=1e-5)
    x1, y1 = curve(t - h)
    x2, y2 = curve(t + h)
    dx = (x2 - x1) / (2h)
    dy = (y2 - y1) / (2h)
    return (dx, dy) ./ norm((dx, dy))
end


function generate_osteocyte_points(N::Int, R::Vector{Float64}, distribution::String; a::Float64 = 1.0, b::Float64 = 1.0)
    curves = [(t,a,b) -> (r * a * cos(2π * t), r * b * sin(2π * t)) for r in R]
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
                x, y = curve(t,a,b)
                push!(points, (x, y))
                push!(distances, norm([x,y]))
                # calculate derivatives
                tangent = numerical_derivative(t -> curve(t,a,b), t; h=1e-5)
                push!(tangents, tangent)
                # Inward normal
                push!(normals, (-tangent[2], tangent[1]))  
            end
        end
    elseif distribution == "Random"
        for _ in 1:N
            curve_idx = rand(1:length(curves))
            curve = curves[curve_idx]
            r = R[curve_idx]
            t = rand()
            x, y = curve(t,a,b)
            push!(points, (x, y))
            push!(distances, norm([x,y]))
             # calculate derivatives
            tangent = numerical_derivative(t -> curve(t,a,b), t; h=1e-5)
            push!(tangents, tangent)
            # Inward normal
            push!(normals, (-tangent[2], tangent[1]))  
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
    end
    # returning Osteocytes storted by distance
    return sort(Osteocytes; by = x -> x.dist_from_origin)
end
