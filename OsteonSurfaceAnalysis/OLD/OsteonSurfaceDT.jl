using FileIO
using ImageIO
using Images
using Colors
using DistanceTransforms
using ImageMorphology
using NPZ
using Printf
using GLMakie

# --- parameters ---
downsample = 1

DIMS = (1024 ÷ downsample, 1024 ÷ downsample)
Z_LAYERS = 202

const BLACK = RGB{N0f8}(0, 0, 0)
const RED   = RGB{N0f8}(1, 0, 0)
const GREEN = RGB{N0f8}(0, 1, 0)

# --- allocate masks ---
outer = falses(DIMS[1], DIMS[2], Z_LAYERS)  # (H, W, Z)
inner = trues(DIMS[1], DIMS[2], Z_LAYERS)
outer_DT = zeros(size(outer))
inner_DT = zeros(size(outer))
outer_DT_inv = zeros(size(outer))
inner_DT_inv = zeros(size(outer))

# --- build masks from slices ---
paths = readdir("./OsteonSurfaceAnalysis/image-masks"; join=true)
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

# --- distance transforms ---
# edt on a Bool array returns distances to the nearest false (background) for true elements.
# To mirror Python's edt on mask and on logical_not(mask), compute both.
for ii in 1:Z_LAYERS
    outer_DT[:,:,ii]     = sqrt.( transform(boolean_indicator(outer[:,:,ii]))  )           # distance inside 'outer' to outside
    inner_DT[:,:,ii]     = sqrt.( transform(boolean_indicator(inner[:,:,ii])) )
    outer_DT_inv[:,:,ii] = sqrt.( transform(boolean_indicator(.!outer[:,:,ii])) )      # distance outside 'outer' to inside
    inner_DT_inv[:,:,ii] = sqrt.( transform(boolean_indicator(.!inner[:,:,ii])) )
end

"""
GLMakie.activate!()
f = Figure(size=(500,500))
ax = CairoMakie.Axis3(f[1,1])

H,W = DIMS
x = collect(1:H)
y = collect(1:W)

ii = 1

surface!(ax,x,y,outer_DT[:,:,ii] .- outer_DT_inv[:,:,ii], colormap=:jet, colorrange=(0,100))
surface!(ax,x,y,inner_DT[:,:,ii] .- inner_DT_inv[:,:,ii], colormap=:jet, colorrange=(0,100))

t = 0.2
values = (1-t).*(outer_DT[:,:,ii] .- outer_DT_inv[:,:,ii])-(t).*(inner_DT[:,:,ii] .- inner_DT_inv[:,:,ii])
surface!(ax,x,y,values, colormap=:jet, colorrange=(0,100))

t = 0.3
values = (1-t).*(outer_DT_inv[:,:,ii])-(t).*(inner_DT_inv[:,:,ii])
map!(x -> (x .< 0 ? 50 : x), values, values)
surface!(ax,x,y, values, colormap=:jet, colorrange=(0,100))

tsamples = 20
dt = 1 / tsamples
tvals = collect(0:dt:1.0)

surface!(ax,x,y,(outer_DT[:,:,ii]), colormap=:jet, colorrange = (0,100))



display(f)
"""

# --- surfaces / phi volume ---
tsamples = 30
dt = 1 / tsamples
tvals = collect(0:dt:1.0)  # length = tsamples + 1

phi = zeros(Float32, DIMS[1], DIMS[2], Z_LAYERS, length(tvals));
values = similar(phi[:, :, 1, 1]);

for (ti, t) in enumerate(tvals)  # ti is 1-based
    println("t index = ", ti, " t value = ", t)
    for ii in 1:Z_LAYERS
        if ti == 1
            phi[:, :, ii, ti] .= outer_DT[:,:,ii] .+ outer_DT_inv[:,:,ii]
        elseif ti == length(tvals)
            phi[:, :, ii, ti] .= inner_DT_inv[:,:,ii] .+ inner_DT[:,:,ii]
        else
            #values .= (1-t).*(outer_DT_inv[:,:,ii])-(t).*(inner_DT_inv[:,:,ii])
            #phi[:, :, ii, ti] .= map!(x -> (x < 0 ? 50 : x), values, values)
            phi[:,:,ii,ti] .= (1-t).*(outer_DT[:,:,ii] .- outer_DT_inv[:,:,ii])-(t).*(inner_DT[:,:,ii] .- inner_DT_inv[:,:,ii])
        end
    end
end

# getting zero contours
import Contour as CTR

H,W = size(phi[:,:,1,1])
x = collect(1:H)
y = collect(1:W)

M = zeros(H,W)

for (ti, t) in enumerate(tvals)
    println("t = ", t)
    vals = findall( x -> (x < 1), phi[:,:,100,ti])
    M[vals] .= 1
end

cset  = CTR.contours(x, y, M, [0])
lvls  = CTR.levels(cset)
lines = CTR.lines(lvls[1])

function Ω(x, y)
    A = 0.0
    n = length(x)
    for ii in 1:n
        j = ii == n ? 1 : ii + 1
        A += x[ii]*y[j] - y[ii]*x[j]
    end
    return abs(A) / 2
end

function find_center(lines)
    # Helper function
    function calc_shape_centroid(x, y)
        A = Ω(x,y)

        x_centroid = 1 / (6 * A) * sum((x[1:end-1] + x[2:end]) .* (x[1:end-1] .* y[2:end] - x[2:end] .* y[1:end-1]))
        y_centroid = 1 / (6 * A) * sum((y[1:end-1] + y[2:end]) .* (x[1:end-1] .* y[2:end] - x[2:end] .* y[1:end-1]))

        return [x_centroid, y_centroid]
    end
    
    A = Inf
    idx = 1
    for (ii,line) in enumerate(lines)
        pts = CTR.coordinates(line)
        x  = pts[1]
        y  = pts[2]
        A_line = Ω(x,y)
        if A_line < A
            A = A_line
            idx = ii
        end
    end

    pts = CTR.coordinates(lines[idx])
    x  = pts[1]
    y  = pts[2]
    return calc_shape_centroid(x, y)
end

function sort_lines_by_radius(lines, center)
    # center: Tuple or SVector (cx, cy)
    radii = [maximum(sqrt.((CTR.coordinates(line)[1] .- center[1]).^2 .+ (CTR.coordinates(line)[2] .- center[2]).^2)) for line in lines]
    sorted_indices = reverse!(sortperm(radii))[1:2:end]
    return lines[sorted_indices]#, radii[sorted_indices]
end

lines_sorted = sort_lines_by_radius(lines, find_center(lines));

GLMakie.activate!()
f = Figure(size=(500,500))
ax = CairoMakie.Axis(f[1,1])
for line in lines_sorted
    pts = CTR.coordinates(line)
    xs  = pts[1]
    ys  = pts[2]
    lines!(ax, xs, ys, linewidth=3, color=:blue)
end

display(f)


function calc_tissue_formation_Rate(interfaces, Δt)
    n = length(interfaces)
    areas = [Ω(CTR.coordinates(interfaces[i])[1], CTR.coordinates(interfaces[i])[2]) for i in 1:n]
    return - diff(areas) ./ Δt
end

calc_tissue_formation_Rate(lines_sorted, dt)

f = Figure(size=(500,500))
ax = CairoMakie.Axis(f[1,1])
lines!(ax,tvals[2:end-2], calc_tissue_formation_Rate(lines_sorted, dt), linewidth=3)
display(f)


function extract_sorted_contours(phi, tvals)

    H, W, Z, _ = size(phi)
    x = collect(1:H)
    y = collect(1:W)
    all_lines_sorted = Vector{Any}(undef, Z)

    for i in 1:Z
        M = zeros(H, W)
        for (ti, t) in enumerate(tvals)
            vals = findall(x -> (x < 1), phi[:, :, i, ti])
            M[vals] .= 1
        end
        cset  = CTR.contours(x, y, M, [0])
        lvls  = CTR.levels(cset)
        lines = CTR.lines(lvls[1])

        center = find_center(lines)
        lines_sorted = sort_lines_by_radius(lines, center)
        all_lines_sorted[i] = lines_sorted
    end

    return all_lines_sorted
end

all_contours = extract_sorted_contours(phi, tvals);

GLMakie.activate!()
f = Figure(size=(500,500))
ax = CairoMakie.Axis(f[1,1])
for line in all_contours[1]
    pts = CTR.coordinates(line)
    xs  = pts[1]
    ys  = pts[2]
    lines!(ax, xs, ys, linewidth=3, color=:blue)
end


display(f)

function estimate_Ot_embed_time(all_contours, Ot_pos, tvals)
    Z_layer = 0;
    # find the closest Z_layer to the z-coordinate of osteocyte
    

    # find the closest contour index to (x,y) position of osteocyte


    return tvals[ii]
end

Ot_pos = [150; 200.5; ]