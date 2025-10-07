# --- Pkg setup (run once) ---
# ] add FileIO ImageIO Images ColorTypes ColorVectorSpace Contour
# ] add GeometryBasics MeshIO

using FileIO, ImageIO, Images
using ColorTypes, ColorVectorSpace
using Contour
using GeometryBasics
using MeshIO
using Statistics
using Interpolations
using GLMakie
using Dierckx

# ---------- helpers ----------
function color_mask(img::AbstractArray{<:Colorant}, which::Symbol)
    RGBimg = convert.(RGB{N0f8}, img)
    R = Float64.(channelview(RGBimg)[1, :, :])
    G = Float64.(channelview(RGBimg)[2, :, :])
    B = Float64.(channelview(RGBimg)[3, :, :])
    which === :red   && return (R .> 0.7) .& (G .< 0.3) .& (B .< 0.3)
    which === :green && return (G .> 0.7) .& (R .< 0.3) .& (B .< 0.3)
    error("which must be :red or :green")
end

"""
Return an ordered boundary polyline for a binary mask as an N×2 matrix (x,y)
in Cartesian coords (origin at bottom-left). Uses Contour.contours(x,y,z,levels).
"""
function mask_boundary_xy(mask::AbstractMatrix{Bool})
    H, W = size(mask)

    # Contour expects z to be sized (length(x), length(y)) = (W, H)
    z = transpose(Float64.(mask))
    x = collect(1:W)          # columns
    y = collect(1:H)          # rows (top→bottom)

    cset  = Contour.contours(x, y, z, [0.5])
    lvls  = Contour.levels(cset)
    isempty(lvls) && return zeros(0,2)

    # Find the longest polyline across all levels (usually there’s just one)
    # there should only be 1 level
    best_line = nothing
    if size(lvls[1].lines,1) > 1
        n_max = 0 
        for ln in Contour.lines(lvls[1])
            n = length(Contour.coordinates(ln)[1])
            if n > n_max
                best_line = ln
                n_max = n
            end
        end
    else
        best_line = Contour.lines(lvls[1])[1]
    end

    #best_line = nothing
    #best_len  = 0
    #for lvl in lvls
    #    for ln in Contour.lines(lvl)
    #        n = length(Contour.coordinates(ln))
    #        if n > best_len
     #           best_len = n
     #           best_line = ln
     #       end
     #   end
    #end
    best_line === nothing && return zeros(0,2)

    pts = Contour.coordinates(best_line)  # vector of points with fields x,y
    xs  = pts[1]
    ys  = pts[2]   # flip y to Cartesian

    # drop duplicate closing point if present
    #if length(xs) > 1 && xs[1] ≈ xs[end] && ys[1] ≈ ys[end]
    #    xs = xs[1:end-1]; ys = ys[1:end-1]
    #end
    hcat(xs, ys)
end

"Resample a closed ring to exactly M points along arclength."
function resample_ring(xy::AbstractMatrix{<:Real}, M::Int)
    @assert size(xy,1) ≥ 3
    P = vcat(xy, xy[1, :]')                           # close
    d = sqrt.(sum(diff(P; dims=1).^2; dims=2))[:]
    s = cumsum(vcat(0.0, d)); total = s[end]
    st = range(0, total, length=M+1)[1:end-1]
    xt = similar(st); yt = similar(st)
    i = 1
    @inbounds for k in eachindex(st)
        t = st[k]
        while s[i+1] < t && i < length(s)-1; i += 1; end
        α = (t - s[i]) / (s[i+1] - s[i] + eps())
        xt[k] = (1-α)*P[i,1] + α*P[i+1,1]
        yt[k] = (1-α)*P[i,2] + α*P[i+1,2]
    end
    hcat(xt, yt)
end

"Ensure points are in a single closed ring ordered by angle; resample to M points."
function ring_resample_angle(xyz::AbstractMatrix{<:Real}, M::Int)
    @assert size(xyz,2) == 3
    x, y = xyz[:,1], xyz[:,2]
    # center for angle param
    cx, cy = mean(x), mean(y)
    θ = atan.(y .- cy, x .- cx)
    # sort by θ (unwrap makes it monotonic)
    p = sortperm(θ)
    θs = θ[p]; xs = x[p]; ys = y[p]
    # enforce 0..2π and uniqueness
    #θu, iu = unique(θs, returninds=true)  # drop duplicates if any
    θu = θs .- minimum(θs)
    θu = 2π .* (θu ./ maximum(θu))        # normalize to [0, 2π]

    # close the ring
    θc = vcat(θu, 2π)
    xc = vcat(xs, xs[1])
    yc = vcat(ys, ys[1])

    # target uniform angles
    θt = range(0, 2π; length=M+1)[1:end-1]

    # 1D linear interpolation along θ
    itpx = interpolate((θc,), xc, Gridded(Linear()))
    itpy = interpolate((θc,), yc, Gridded(Linear()))
    xrt = itpx.(θt)
    yrt = itpy.(θt)
    hcat(xrt, yrt)
end

"Rotate second ring so its start index best matches the first (minimize seam twist)."
function align_ring!(ring::AbstractMatrix, ref::AbstractMatrix)
    d = sum((ref[1:1,:] .- ring).^2; dims=2)[:]
    k = argmin(d)
    k == 1 && return ring
    vcat(ring[k:end, :], ring[1:k-1, :])
end

function Ω(p)
    A = 0
    for ii in axes(p,2)
        if ii == size(p,2)
            A += (p[1,ii]*p[2,1] -  p[2,ii]*p[1,1])
        else
            A += (p[1,ii]*p[2,ii+1] -  p[2,ii]*p[1,ii+1])
        end
    end
    return abs(A)/2;
end

function ωκ(rᵢ₋₁, rᵢ, rᵢ₊₁)
    triVector = [rᵢ₋₁' rᵢ' rᵢ₊₁']
    A = zeros(size(triVector,1))
    for ii in axes(triVector,1)
        A[ii] = Ω(reshape(triVector[ii,:],(2,3)))
    end
    return A
end

# Approximating curvature of 2D lines
function κ(rᵢ₋₁, rᵢ, rᵢ₊₁)

    A = ωκ(rᵢ₋₁,rᵢ,rᵢ₊₁)
    l1 = .√(sum((rᵢ- rᵢ₋₁).^2,dims=1))
    l2 = .√(sum((rᵢ₊₁- rᵢ).^2,dims=1))
    l3 = .√(sum((rᵢ₊₁- rᵢ₋₁).^2,dims=1))

    return (4*A)./(l1.*l2.*l3)'
end

# ---------- build interpolants ----------
"""
build_surface_interpolants(contours; M=256)

Input:
  contours :: Vector{Matrix{Float64}}  # each N×3 (x,y,z), one ring per slice

Output:
  (itpX, itpY, zgrid, thetagrid)

Usage:
  xs = itpX(z0, θ)   # z0 scalar, θ scalar
  xs = itpX.(z0, θs) # broadcast over θs
"""
function build_surface_interpolants(contours::Vector{Matrix{Float64}}; M::Int=256)
    @assert !isempty(contours)
    # sort slices by z (assumes each contour row shares same z)
    zs = [mean(C[:,3]) for C in contours]
    order = sortperm(zs)
    contours = contours[order]
    zgrid = sort(zs)

    # resample each to M angles and align seams
    rings = Matrix{Float64}[]
    for (i, C) in enumerate(contours)
        r = ring_resample_angle(C, M)  # M×2 (x,y)
        if i > 1
            r = align_ring!(r, rings[end])
        end
        push!(rings, r)
    end
    Nz = length(rings)

    # build data arrays X(z,θ), Y(z,θ) with shape (Nz, M)
    X = Array{Float64}(undef, Nz, M)
    Y = Array{Float64}(undef, Nz, M)
    for s in 1:Nz
        X[s, :] = rings[s][:,1]
        Y[s, :] = rings[s][:,2]
    end

    # angular axis (uniform)
    thetagrid = range(0, 2π; length=M+1)[1:end-1]

    # make θ periodic by appending a wrap column at 2π (duplicate of column 1)
    Xp = hcat(X, X[:,1])
    Yp = hcat(Y, Y[:,1])
    thetap = vcat(collect(thetagrid), 2π)

    # Build gridded linear interpolants over (z, θ)
    itpX = interpolate((zgrid, thetap), Xp, Gridded(Linear()))
    itpY = interpolate((zgrid, thetap), Yp, Gridded(Linear()))
    #itpX = interpolate((zgrid, thetap), Xp, BSpline(Constant()))
    #itpY = interpolate((zgrid, thetap), Yp, BSpline(Constant()))

    return [itpX, itpY, zgrid, thetagrid]
end


# ---------- main ----------
"""
build_surfaces(paths; Δz=0.1, M=200) -> (outer_mesh, inner_mesh)
`paths`: image files ordered by slice (z ascending). `M`: samples per ring.
"""
function build_surfaces(paths::Vector{String}; Δz::Float64=0.1, M::Int=200)
    outer_rings = Matrix{Float64}[]
    inner_rings = Matrix{Float64}[]
    #curvature_outer = Matrix{Float64}[]
    #curvature_inner = Matrix{Float64}[]
    for (i,p) in enumerate(paths)
        img = load(p)
        r_xy = mask_boundary_xy(color_mask(img, :red))
        g_xy = mask_boundary_xy(color_mask(img, :green))
        isempty(r_xy) && error("No red boundary in $p")
        isempty(g_xy) && error("No green boundary in $p")

        push!(outer_rings, hcat([resample_ring(r_xy, M),ones(size(resample_ring(r_xy, M)[:,1])).*((i-1)*Δz)]...))
        push!(inner_rings, hcat([resample_ring(g_xy, M),ones(size(resample_ring(g_xy, M)[:,1])).*((i-1)*Δz)]...))
    end
    
    outer = build_surface_interpolants(outer_rings)
    inner = build_surface_interpolants(inner_rings)

    # calculating curvature of each z ring 
    curvature_outer = similar(outer[1].coefs)
    curvature_inner = similar(inner[1].coefs)
    for ii in 1:size(inner[1].coefs,1)
        x_inner = inner[1].coefs[ii,:]
        y_inner = inner[2].coefs[ii,:]
        κ_inner = κ(circshift(hcat(x_inner,y_inner)'[:,1:end-1], (0,1)), hcat(x_inner,y_inner)'[:,1:end-1], circshift(hcat(x_inner,y_inner)'[:,1:end-1], (0,-1)))
        x_outer = outer[1].coefs[ii,:]
        y_outer = outer[2].coefs[ii,:]
        κ_outer = κ(circshift(hcat(x_outer,y_outer)'[:,1:end-1], (0,1)), hcat(x_outer,y_outer)'[:,1:end-1], circshift(hcat(x_outer,y_outer)'[:,1:end-1], (0,-1)))

        curvature_outer[ii,:] = [κ_outer;κ_outer[1]]
        curvature_inner[ii,:] = [κ_inner;κ_inner[1]]
    end

    return outer, inner, outer_rings, inner_rings, curvature_outer, curvature_inner
end



# Example:
#paths = sort(readdir("path/to/masks"; join=true))
#paths = ["./image-masks/FM40-10-mask0000.png", "./image-masks/FM40-10-mask0200.png","./image-masks/FM40-10-mask0201.png"]#, "./image-masks/FM40-10-mask0002.png", "./image-masks/FM40-10-mask0003.png", "./image-masks/FM40-10-mask0004.png","./image-masks/FM40-10-mask0005.png","./image-masks/FM40-10-mask0006.png","./image-masks/FM40-10-mask0007.png", "./image-masks/FM40-10-mask0008.png", "./image-masks/FM40-10-mask0009.png", "./image-masks/FM40-10-mask0010.png"]

paths = readdir("./OsteonSurfaceAnalysis/image-masks"; join=true)
outer, inner, outer_rings, inner_rings, curvature_outer, curvature_inner = build_surfaces(paths; Δz=0.0625, M=400)

x_outer = outer[1].coefs
y_outer = outer[2].coefs
z_outer = repeat(outer[3], outer=(1, size(outer[1],2)))


x_inner = inner[1].coefs
y_inner = inner[2].coefs
z_inner = repeat(inner[3], outer=(1, size(inner[1],2)))


GLMakie.activate!()
fig = GLMakie.Figure(size = (800, 600))
ax = GLMakie.Axis3(fig[1, 1], perspectiveness = 0.75, xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false, xlabelvisible = false, ylabelvisible = false, zlabelvisible = false)
GLMakie.surface!(ax, x_outer, y_outer, z_outer, color = curvature_outer, colormap = :plasma, colorrange = (0, 0.1))
GLMakie.surface!(ax, x_inner, y_inner, z_inner, color = curvature_inner, colormap = :dense, colorrange = (0, 1.0))
#GLMakie.Colorbar(fig[1,2], colormap=:jet, colorrange=(-0.02, 0.7), label = "z slice curvature", width = 20)

GLMakie.lines!(ax, outer_rings[1][:,1], outer_rings[1][:,2], outer_rings[1][:,3], color = :red)
GLMakie.lines!(ax, inner_rings[1][:,1], inner_rings[1][:,2], inner_rings[1][:,3], color = :blue)


start_angle = π / 4
n_frames = 200
ax.viewmode = :fit # Prevent axis from resizing during animation
record(fig, "mask_to_3D_Osteons.mp4", 1:n_frames) do frame
    ax.azimuth[] = start_angle + 2pi * frame / n_frames
end