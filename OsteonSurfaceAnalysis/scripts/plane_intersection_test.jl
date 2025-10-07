# ----- tuple helpers -----
dot3(a,b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
minus(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
plus(a,b)  = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
scale(a,s) = (a[1]*s, a[2]*s, a[3]*s)
cross3(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
norm3(a) = sqrt(dot3(a,a))
normalize3(v) = (s = norm3(v); s==0 ? (0.0,0.0,0.0) : scale(v, 1/s))

struct Plane{T}
    p0::NTuple{3,T}   # a point on the plane
    n::NTuple{3,T}    # (not necessarily unit) normal
end

# Rodrigues rotation of a vector v around unit axis û by angle θ
function rotate_about_axis(v::NTuple{3,Float64}, û::NTuple{3,Float64}, θ::Float64)
    c, s = cos(θ), sin(θ)
    term1 = scale(v, c)
    term2 = scale(cross3(û, v), s)
    term3 = scale(û, dot3(û, v) * (1 - c))
    plus(plus(term1, term2), term3)
end


"""
    plane_through_centers(top, bot; θ=0.0, ref=(0.0,0.0,1.0))

Return a plane `Plane(p0, n)` that:
- contains the line through `top` and `bot`
- is obtained by rotating an initial plane around that line by angle `θ` (radians)

`ref` picks the initial orientation before rotation (any vector not parallel to the line).
"""
function plane_through_centers(top::NTuple{3,Float64}, bot::NTuple{3,Float64};
                               θ::Float64=0.0, ref::NTuple{3,Float64}=(0.0,0.0,1.0))
    # Axis of rotation = line through the centers
    axis = minus(top, bot)
    û = normalize3(axis)
    if norm3(û) == 0
        error("Top and bottom centers coincide; axis undefined.")
    end

    # Build a vector inside the plane that is orthogonal to the axis:
    # take ref, remove its component along the axis
    v0 = minus(ref, scale(û, dot3(ref, û)))
    if norm3(v0) == 0
        # ref was parallel to the axis; pick another
        ref2 = abs(û[3]) < 0.9 ? (0.0,0.0,1.0) : (1.0,0.0,0.0)
        v0 = minus(ref2, scale(û, dot3(ref2, û)))
    end
    v0 = normalize3(v0)

    # Initial plane normal is n0 = axis × v0  (so plane contains the axis)
    n0 = cross3(û, v0)

    # Rotate the plane normal around the axis by θ
    nθ = rotate_about_axis(n0, û, θ)

    # Plane through any point on the line (use top)
    return Plane{Float64}(top, nθ)
end

function intersect_segment_with_plane(p::NTuple{3,Float64},
                                      q::NTuple{3,Float64},
                                      pl::Plane{Float64}; eps=1e-12)
    pq = minus(q, p)
    denom = dot3(pl.n, pq)
    num   = -dot3(pl.n, minus(p, pl.p0))
    if abs(denom) < eps
        return false, (0.0,0.0,0.0)   # parallel/coplanar case ignored
    end
    t = num/denom
    if t < -eps || t > 1+eps
        return false, (0.0,0.0,0.0)
    end
    return true, plus(p, scale(pq, clamp(t, 0.0, 1.0)))
end

function intersect_polylines_with_plane(polys::Vector{Vector{NTuple{3,Float64}}},
                                        pl::Plane{Float64}; closed::Bool=false)
    hits = NTuple{3,Float64}[]
    for poly in polys
        n = length(poly); n < 2 && continue
        segs = closed ? [(i, i % n + 1) for i in 1:n] : [(i, i+1) for i in 1:n-1]
        for (i,j) in segs
            hit, X = intersect_segment_with_plane(poly[i], poly[j], pl)
            hit && push!(hits, X)
        end
    end
    hits
end

# From your code:
# ϕ_top = ϕ[:, :, zidx_top, tidx]
# cset_top = CTR.contours(x, y, ϕ_top, [0.0])

function contours_to_3d_polylines(cset, zval)
    polys = Vector{Vector{NTuple{3,Float64}}}()
    for lvl in CTR.levels(cset)
        for ln in CTR.lines(lvl)
            xs, ys = CTR.coordinates(ln)
            push!(polys, [(xs[i], ys[i], zval) for i in eachindex(xs)])
        end
    end
    polys
end

# centers you already computed from inner-most contours:
top_center    = (cx_top, cy_top, Float64(z_top))       # Float64 tuples
bottom_center = (cx_bot, cy_bot, Float64(z_bot))

# Build plane that contains the line top↔bottom; rotate by θ radians about that line
θ = deg2rad(90)  # example rotation
pl = plane_through_centers(top_center, bottom_center; θ)

# Build 3D polylines for all relevant contours (from top & bottom slices, or multiple z’s)
polys3d = vcat(contours_to_3d_polylines(cset_top, z_top),
               contours_to_3d_polylines(cset_bot, z_bot))

# Intersections between those 3D polylines and your rotated plane
hits3d = intersect_polylines_with_plane(polys3d, pl; closed=true)

# `hits3d` are the 3D points where the contour lines meet your rotated plane.
scatter!(a1,hits3d[1][1], hits3d[1][2], hits3d[1][3], markersize=25)
scatter!(a1,hits3d[2][1], hits3d[2][2], hits3d[2][3], markersize=25)
scatter!(a1,hits3d[3][1], hits3d[3][2], hits3d[3][3], markersize=25)
scatter!(a1,hits3d[4][1], hits3d[4][2], hits3d[4][3], markersize=25)