# Example usage
Tdelay = 0.5
path_to_processed_img = readdir("./OsteonSurfaceAnalysis/DATA/FM40-1-R2/Processed_Images"; join=true)[end] 
path_to_HCa = readdir("./OsteonSurfaceAnalysis/DATA/FM40-1-R2/HCa"; join = true)[end]
path_to_output = "./OsteonSurfaceAnalysis/scripts/test2.png"

make_tdelay_mask(Tdelay,path_to_processed_img,path_to_HCa,path_to_output)

path_to_bottom = readdir("./OsteonSurfaceAnalysis/DATA/FM40-1-R2/Processed_Images"; join=true)[1]
path_to_top = path_to_output

paths = [path_to_bottom, path_to_top]

# --- parameters ---
downsample = 1

DIMS = (1024 ÷ downsample, 1024 ÷ downsample)
Z_LAYERS = 2

# --- allocate masks ---
outer = falses(DIMS[1], DIMS[2], Z_LAYERS);  # (H, W, Z)
inner = trues(DIMS[1], DIMS[2], Z_LAYERS);
outer_DT = zeros(size(outer));
inner_DT = zeros(size(outer));
outer_DT_inv = zeros(size(outer));
inner_DT_inv = zeros(size(outer));

for z0 in 0:Z_LAYERS-1
    fn = paths[z0+1]
    img = load(fn)
     # flip vertically
    #img = reverse(img, dims=1)
    img = reverse(img, dims=2)

    # downsample
    img = img[1:downsample:end, 1:downsample:end]

    # color-based masks
    m_green = (img .== GREEN)
    m_red   = (img .== RED)

    outer[:, :, z0+1] .= m_green .| m_red
    inner[:, :, z0+1] .= .!m_green
end

outer_dt_S = zeros(size(outer,1), size(outer,2),Z_LAYERS);
inner_dt_S = zeros(size(outer,1), size(outer,2),Z_LAYERS);

for z0 in 1:Z_LAYERS
    outer_dt_S[:,:,z0] = edt_S(outer[:,:,z0]);
    inner_dt_S[:,:,z0] = edt_S(inner[:,:,z0]);
end

tsamples = 5
dt = 1 / tsamples
tvals = collect(0:dt:1.0)  # length = tsamples + 1

ϕ_func = (t,S_DTʰ,S_DTᶜ) -> (1-t) .* S_DTʰ - (t) .* S_DTᶜ

ϕ = zeros(Float32, DIMS[1], DIMS[2], Z_LAYERS, length(tvals));

# 2D level set function
for (ti, t) in enumerate(tvals)  # ti is 1-based
    println("t index = ", ti, " t value = ", t)
    for z in 1:Z_LAYERS 
        ϕ[:,:,z,ti] .= ϕ_func(t, outer_dt_S[:,:,z], inner_dt_S[:,:,z]) 
    end
end


H,W,D = size(ϕ[:,:,:,1])
x = collect(1:H)
y = collect(1:W)
Δz = 50.0;

set_theme!(theme_black(), fontsize = 36)

GLMakie.activate!()
f = Figure(size = (800, 800))

a1 = Axis3(f[1, 1], title = "ϕ contours")

for (t,ti) in enumerate(tvals)
    ϕ_bottom = ϕ[:,:,1,t]
    cset_bot = CTR.contours(x,y,ϕ_bottom, [0])
    line_bot = first(CTR.lines(first(CTR.levels(cset_bot))))
    x_bot, y_bot = CTR.coordinates(line_bot)

    ϕ_top = ϕ[:,:,2,t]
    cset_top = CTR.contours(x,y,ϕ_top, [0])
    line_top = first(CTR.lines(first(CTR.levels(cset_top))))
    x_top, y_top = CTR.coordinates(line_top)

    lines!(a1, x_bot, y_bot, zeros(size(x_bot)), linewidth=3, color=:red)
    lines!(a1, x_top, y_top, Δz.*ones(size(x_top)), linewidth=3, color=:blue)
end

function Ω(x, y)
    A = 0.0
    n = length(x)
    for ii in 1:n
        j = ii == n ? 1 : ii + 1
        A += x[ii]*y[j] - y[ii]*x[j]
    end
    return abs(A) / 2
end

function find_center(x, y)
    # Ensure the polygon is closed
    if x[1] != x[end] || y[1] != y[end]
        x = vcat(x, x[1])
        y = vcat(y, y[1])
    end

    A = Ω(x, y)

    x_centroid = 1 / (6 * A) * sum((x[1:end-1] + x[2:end]) .* (x[1:end-1] .* y[2:end] - x[2:end] .* y[1:end-1]))
    y_centroid = 1 / (6 * A) * sum((y[1:end-1] + y[2:end]) .* (x[1:end-1] .* y[2:end] - x[2:end] .* y[1:end-1]))

    return [x_centroid, y_centroid]
end

# Functions to generate z-cutting planes and calculate which points from the contours intersect with the plane
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

# calculating centers

ϕ_bottom = ϕ[:,:,1,6]
cset_bot = CTR.contours(x,y,ϕ_bottom, [0])
line_bot = first(CTR.lines(first(CTR.levels(cset_bot))))
x_bot, y_bot = CTR.coordinates(line_bot)
center_bot = find_center(reverse(x_bot), reverse(y_bot))

ϕ_top = ϕ[:,:,2,6]
cset_top = CTR.contours(x,y,ϕ_top, [0])
line_top = first(CTR.lines(first(CTR.levels(cset_top))))
x_top, y_top = CTR.coordinates(line_top)
center_top = find_center(x_top, y_top)

scatter!(a1, center_bot[1], center_bot[2], 0, markersize=30)
scatter!(a1, center_top[1], center_top[2], Δz, markersize=30)

θs = collect(0.0:2π/5:2π)
top_center = (center_top[1], center_top[2], Δz)
bottom_center = (center_bot[1], center_bot[2], 0.0)

intersecting_points_per_contour = []

for (ti,t) in enumerate(tvals)
    println(ti)
    ϕ_bottom = ϕ[:,:,1,ti]
    ϕ_top = ϕ[:,:,2,ti]
    cset_top = CTR.contours(x,y,ϕ_top, [0])
    cset_bot = CTR.contours(x,y,ϕ_bottom, [0])
    intersecting_points_per_theta = []
    for θ in θs
        # Build plane that contains the line top↔bottom; rotate by θ radians about that line
        pl = plane_through_centers(top_center, bottom_center; θ)
        # Build 3D polylines for all relevant contours (from top & bottom slices, or multiple z’s)
        polys3d = vcat(contours_to_3d_polylines(cset_top, z_top),
               contours_to_3d_polylines(cset_bot, z_bot))
        hits3d = intersect_polylines_with_plane(polys3d, pl; closed=true)
        push!(intersecting_points_per_theta, hits3d)
    end
    push!(intersecting_points_per_contour, intersecting_points_per_theta)
end

for pts_in_contour in intersecting_points_per_contour
    for points_per_angle in pts_in_contour
        scatter!(a1,points_per_angle[1][1], points_per_angle[1][2], points_per_angle[1][3], markersize=25, color=:green)
        scatter!(a1,points_per_angle[2][1], points_per_angle[2][2], points_per_angle[2][3], markersize=25, color=:yellow)
        scatter!(a1,points_per_angle[3][1], points_per_angle[3][2], points_per_angle[3][3], markersize=25, color=:purple)
        scatter!(a1,points_per_angle[4][1], points_per_angle[4][2], points_per_angle[4][3], markersize=25, color=:cyan)
    end
end

