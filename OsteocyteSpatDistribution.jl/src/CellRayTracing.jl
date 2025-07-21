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

function generate_rays_from_cell(cell::Ellipse, max_ray_length=10, num_rays=100)
    θs = range(0, 2π, length=num_rays+1)[1:end-1]  # Uniform angles, exclude endpoint
    rays = []
    for θ in θs
        dir = [cos(θ), sin(θ)]
        ray_start = cell.center
        ray_end = cell.center + max_ray_length * dir
        push!(rays, [ray_start, ray_end])
    end
    return rays
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

function calc_ray_lengths(intersecting_rays::Vector{Any})
    ray_lengths = []
    for rays_from_cell in intersecting_rays
        for ray in rays_from_cell
            length = norm(ray[2] - ray[1])
            push!(ray_lengths, length)
        end
    end
    return ray_lengths
end

function acos2(x, y)
    # Returns the signed angle between vectors x and y, in [-π, π]
    angle = acos(clamp(dot(x, y) / (norm(x) * norm(y)), -1.0, 1.0))
    sign = signbit(x[1]*y[2] - x[2]*y[1]) ? 1 : -1
    return angle * sign
end

function calc_ray_angle(cell::Ellipse, rays_from_cell::Vector{Any})
    angles = []
    # Calculate vector of cell long axis
    a,b = cell.axes
    R = [cos(cell.angle) -sin(cell.angle); sin(cell.angle) cos(cell.angle)]
    long_axis_point = [R[1] * a * cos(0), R[2] * b * sin(0)] 
    cell_vec = long_axis_point - cell.center
    for ray in rays_from_cell
        ray_start, ray_end = ray
        ray_vec = ray_end - ray_start
        # Calculate angle with respect to the cell's orientation
        angle = acos2(ray_vec, cell_vec)
        push!(angles, angle)
    end
    return angles
end

function calc_ray_angles(cells::Vector{Ellipse}, intersecting_rays::Vector{Any})
    all_angles = []
    for (i,cell) in enumerate(cells)
        angles = calc_ray_angle(cell, intersecting_rays[i])
        push!(all_angles, angles)
    end
    return all_angles
end

function separate_forward_backward_rays(cells::Vector{Ellipse}, intersecting_rays::Vector{Any})
    forward_rays = []
    backward_rays = []
    ray_angles = calc_ray_angles(cells, intersecting_rays)
    for (i,cell_rays) in enumerate(ray_angles)
        forward_rays_from_cell = []
        backward_rays_from_cell = []
        for (j, ray) in enumerate(cell_rays)
            if ray < 0
                push!(backward_rays_from_cell, intersecting_rays[i][j])
            else
                push!(forward_rays_from_cell, intersecting_rays[i][j])
            end
        end
        push!(forward_rays, forward_rays_from_cell)
        push!(backward_rays, backward_rays_from_cell)
    end
    return forward_rays, backward_rays
end

using OrderedCollections

function count_ray_intersections(cells::Vector{Ellipse}, intersecting_rays::Vector{Any})
    intersection_counts = []
    for (i, cell) in enumerate(cells)
        counts = OrderedDict{Int, Int}()
        # Initialize all other indices with 0
        for j in 1:length(cells)
            if j != i
                counts[j] = 0
            end
        end
        rays = intersecting_rays[i]
        for ray in rays
            # Use destructuring to assign and use ray_start and ray_end
            ray_start = ray[1]
            ray_end = ray[2]
            for (j, other_cell) in enumerate(cells)
                if j != i
                    other_pts = generate_cell_from_Ellipse(other_cell)
                    for k in 1:length(other_pts)-1
                        seg_start = other_pts[k]
                        seg_end = other_pts[k+1]
                        seg_vec = seg_end - seg_start
                        pt_vec = ray_end - seg_start
                        seg_len2 = dot(seg_vec, seg_vec)
                        if seg_len2 == 0
                            continue
                        end
                        proj = dot(pt_vec, seg_vec) / seg_len2
                        if 0.0 ≤ proj ≤ 1.0
                            closest_pt = seg_start + proj * seg_vec
                            # Check if the ray_end is very close to the segment (boundary of other cell)
                            if norm(closest_pt - ray_end) < 1e-3
                                counts[j] += 1
                                break
                            end
                        end
                    end
                end
            end
        end
        # Store as OrderedDict for each cell
        push!(intersection_counts, counts)
    end
    return intersection_counts
end