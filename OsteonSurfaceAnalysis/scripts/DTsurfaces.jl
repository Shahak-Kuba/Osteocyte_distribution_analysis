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
Z_LAYERS = 183

const BLACK = RGB{N0f8}(0, 0, 0)
const RED   = RGB{N0f8}(1, 0, 0)
const GREEN = RGB{N0f8}(0, 1, 0)

# --- allocate masks ---
outer = falses(DIMS[1], DIMS[2], Z_LAYERS);  # (H, W, Z)
inner = trues(DIMS[1], DIMS[2], Z_LAYERS);
outer_DT = zeros(size(outer));
inner_DT = zeros(size(outer));
outer_DT_inv = zeros(size(outer));
inner_DT_inv = zeros(size(outer));

# --- build masks from slices ---
paths = readdir("./OsteonSurfaceAnalysis/DATA/FM40-1-R1/Processed_Images"; join=true)
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

function edt(mask::BitArray)
    return sqrt.(DistanceTransforms.transform(boolean_indicator(mask)))
end

function edt_S(mask::BitArray)
    return edt(mask) .- edt(.!mask)
end

outer_dt_S = edt_S(outer);
inner_dt_S = edt_S(inner);

"""
outer_DT  = edt(outer)
inner_DT = edt(inner)
outer_DT_inv = edt(.!outer)      
inner_DT_inv = edt(.!inner)

signed_DT_H = outer_DT .- outer_DT_inv
signed_DT_C = inner_DT .- inner_DT_inv
"""

# --- surfaces / phi volume ---
tsamples = 5
dt = 1 / tsamples
tvals = collect(0:dt:1.0)  # length = tsamples + 1

ϕ_func = (t,S_DTʰ,S_DTᶜ) -> (1-t) .* S_DTʰ - (t) .* S_DTᶜ

ϕ = zeros(Float32, DIMS[1], DIMS[2], Z_LAYERS, length(tvals));

for (ti, t) in enumerate(tvals)  # ti is 1-based
    println("t index = ", ti, " t value = ", t)
    ϕ[:,:,:,ti] .= ϕ_func(t, outer_dt_S, inner_dt_S) #(1-t).*outer_dt_S - (t).*inner_dt_S
    """
    if ti == 1
        phi[:, :, :, ti] .= outer_DT[:,:,ii] .- outer_DT_inv[:,:,ii]
    elseif ti == length(tvals)
        phi[:, :, :, ti] .= inner_DT_inv[:,:,ii] .- inner_DT[:,:,ii]
    else
        phi[:,:,:,ti] .= (1-t).*(outer_DT[:,:,ii] .- outer_DT_inv[:,:,ii])-(t).*(inner_DT[:,:,ii] .- inner_DT_inv[:,:,ii])
    end
    """
end

"""
    kappa = curvature_central(phi, dx, dy, dz; eps=1e-12)

Compute mean curvature κ of the level set φ(x,y,z) on a regular 3-D grid,
using 2nd-order central differences.

κ = ( φx^2 φyy − 2 φx φy φxy + φy^2 φxx
    + φx^2 φzz − 2 φx φz φxz + φz^2 φxx
    + φy^2 φzz − 2 φy φz φyz + φz^2 φyy ) / |∇φ|^3

`phi` must be a 3-D `AbstractArray{<:Real,3}`. The spacing (`dx,dy,dz`)
are the grid steps in x,y,z. Boundaries are set to `NaN` (edit as needed).
"""
function curvature_central(ϕ::AbstractArray{<:Real,3},
                           dx::Real, dy::Real, dz::Real; eps=1e-12)

    nx, ny, nz = size(ϕ)
    kappa = fill!(similar(ϕ, Float64), 0.0)

    inv2dx = 1.0/(2*dx); inv2dy = 1.0/(2*dy); inv2dz = 1.0/(2*dz)
    invdx2 = 1.0/(dx*dx); invdy2 = 1.0/(dy*dy); invdz2 = 1.0/(dz*dz)
    inv4dxdy = 1.0/(4*dx*dy); inv4dxdz = 1.0/(4*dx*dz); inv4dydz = 1.0/(4*dy*dz)

    @inbounds for k in 2:nz-1, j in 2:ny-1, i in 2:nx-1
        ϕc = ϕ[i, j, k]

        # First derivatives (central)
        ϕx = (ϕ[i+1, j,   k  ] - ϕ[i-1, j,   k  ]) * inv2dx
        ϕy = (ϕ[i,   j+1, k  ] - ϕ[i,   j-1, k  ]) * inv2dy
        ϕz = (ϕ[i,   j,   k+1] - ϕ[i,   j,   k-1]) * inv2dz

        # Second derivatives (central)
        ϕxx = (ϕ[i+1, j,   k  ] - 2ϕc + ϕ[i-1, j,   k  ]) * invdx2
        ϕyy = (ϕ[i,   j+1, k  ] - 2ϕc + ϕ[i,   j-1, k  ]) * invdy2
        ϕzz = (ϕ[i,   j,   k+1] - 2ϕc + ϕ[i,   j,   k-1]) * invdz2

        # Mixed derivatives (central)
        ϕxy = (ϕ[i+1, j+1, k] - ϕ[i+1, j-1, k] - ϕ[i-1, j+1, k] + ϕ[i-1, j-1, k]) * inv4dxdy
        ϕxz = (ϕ[i+1, j, k+1] - ϕ[i+1, j, k-1] - ϕ[i-1, j, k+1] + ϕ[i-1, j, k-1]) * inv4dxdz
        ϕyz = (ϕ[i, j+1, k+1] - ϕ[i, j+1, k-1] - ϕ[i, j-1, k+1] + ϕ[i, j-1, k-1]) * inv4dydz

        # |∇ϕ|
        gradmag = sqrt(ϕx*ϕx + ϕy*ϕy + ϕz*ϕz) + eps  # eps avoids divide-by-zero
        denom = gradmag^3

        # Numerator (your formula, grouped by pairs)
        num  =  (ϕx^2)*ϕyy - 2*ϕx*ϕy*ϕxy + (ϕy^2)*ϕxx
        num +=  (ϕx^2)*ϕzz - 2*ϕx*ϕz*ϕxz + (ϕz^2)*ϕxx
        num +=  (ϕy^2)*ϕzz - 2*ϕy*ϕz*ϕyz + (ϕz^2)*ϕyy

        kappa[i, j, k] = num / denom
    end

    return kappa
end

dx = 1; dy = 1; dz = 1;
K = curvature_central(ϕ[:,:,:,5], dx, dy, dz)



H,W,D = size(ϕ[:,:,:,1])
x = collect(1:H)
y = collect(1:W)

set_theme!(theme_black(), fontsize = 36)

GLMakie.activate!()
f = Figure(size = (800, 800))

a1 = Axis3(f[1, 1], title = "ϕ contours")

surface!(a1, x, y, K[:,:,182],colormap=:jet,colorrange=(0.0,0.1))
contour!(a1, 1 .. H, 1 .. W, 1 .. D, ϕ[:,:,:,5], levels = [0], colormap=:jet, colorrange=(0,400))

for ii in 1:length(tvals)
    contour!(a1, 1 .. H, 1 .. W, 1 .. D, ϕ[:,:,:,ii], levels = [0], colormap=:jet, colorrange=(0,400))
end

function calc_t_from_Ot_pos(Ot_pos, DT_H, DT_C)
    x = Ot_pos[1]; y = Ot_pos[2]; z = Ot_pos[3];
    return Float64(DT_H[x,y,z] / (DT_H[x,y,z] - DT_C[x,y,z]))
end

Ot_pos = [400;200;4]

calc_t_from_Ot_pos(Ot_pos,outer_dt_S,inner_dt_S)