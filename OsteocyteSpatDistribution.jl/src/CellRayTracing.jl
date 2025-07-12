
function generate_cell_from_Ellipse(e::Ellipse, num_boundary_points=100)
    θ = range(0, 2π, length=num_boundary_points)
    a, b = e.axes
    # Parametric ellipse
    xs = a * cos.(θ)
    ys = b * sin.(θ)
    # Rotate and translate
    R = [cos(e.angle) -sin(e.angle); sin(e.angle) cos(e.angle)]
    pts = [R * [xs[i], ys[i]] + e.center for i in 1:length(θ)]
    return [[p[1], p[2]] for p in pts]
end

function generate_rays_from_cell(cell::Ellipse, max_ray_length=10)
    cell_bp = generate_cell_from_Ellipse(cell)
    all_rays = []
    for point in cell_bp[1:end-1]
        ray_dir = point - cell.center
        ray = [point, point + (ray_dir / norm(ray_dir)) * max_ray_length]
        push!(all_rays, ray)
    end
    return all_rays
end

function store_intersecting_rays(cells::Vector{Ellipse}, max_distance=50, R=1.0)
    all_intersecting_rays = []
    for cell in cells
        rays = generate_rays_from_cell(cell, max_distance)
        intersecting_rays = []
        for ray in rays
            ray_start, ray_end = ray
            ray_vec = ray_end - ray_start
            closest_t = Inf
            intersection_point = nothing
            for other_cell in cells
                if other_cell != cell
                    other_pts = generate_cell_from_Ellipse(other_cell)
                    for i in 1:length(other_pts)-1
                        seg_start = other_pts[i]
                        seg_end = other_pts[i+1]
                        seg_vec = seg_end - seg_start
                        # Ray-segment intersection
                        function ray_segment_intersect(p, r, q, s)
                            cross_rs = r[1]*s[2] - r[2]*s[1]
                            if abs(cross_rs) < 1e-8
                                return nothing # Parallel
                            end
                            qp = q - p
                            t = (qp[1]*s[2] - qp[2]*s[1]) / cross_rs
                            u = (qp[1]*r[2] - qp[2]*r[1]) / cross_rs
                            if (t ≥ 0) && (u ≥ 0) && (u ≤ 1)
                                return t
                            else
                                return nothing
                            end
                        end
                        t = ray_segment_intersect(ray_start, ray_vec, seg_start, seg_vec)
                        if t !== nothing && t < closest_t && t ≤ 1
                            closest_t = t
                            intersection_point = ray_start + t * ray_vec
                        end
                    end
                end
            end
            if intersection_point !== nothing
                # Check if the ray segment [ray_start, intersection_point] passes within the forbidden disc
                center = [0.0, 0.0]  # Assuming the center of the disc is at the origin
                v = intersection_point - ray_start
                v_norm2 = dot(v, v)
                passes_disc = false
                # Sample points along the segment and check distance to center
                for s in 0:0.05:1
                    pt = ray_start + s * v
                    if norm(pt - center) < R
                        passes_disc = true
                        break
                    end
                end
                if !passes_disc
                    push!(intersecting_rays, [ray_start, intersection_point])
                end
            end
        end
        push!(all_intersecting_rays, intersecting_rays)
    end  
    return all_intersecting_rays
end

intersecting_rays = store_intersecting_rays(cells)

for rays_from_cell in intersecting_rays
    for ray in rays_from_cell
        x = [ray[1][1], ray[2][1]]
        y = [ray[1][2], ray[2][2]]
        lines!(ax, x, y, color=:red, alpha=0.5)
    end
end