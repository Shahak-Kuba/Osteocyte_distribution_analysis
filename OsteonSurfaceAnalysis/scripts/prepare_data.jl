# save as make_mask.jl and run with: julia make_mask.jl
using FileIO, Images, ImageBinarization, ImageMorphology, ImageSegmentation, ImageFiltering, Statistics
using DistanceTransforms
import Contour as CTR

function process_mask_data(path_HCa, path_On, output_path)

    # --- load & prep ---
    hca = load(path_HCa)
    on  = load(path_On)

    # light smoothing to denoise before thresholding
    hca_blur = imfilter(hca, Kernel.gaussian(2))
    on_blur  = imfilter(on,  Kernel.gaussian(2))

    # --- automatic thresholding (Otsu) ---
    HCa = binarize(hca_blur, Otsu())             
    On = binarize(on_blur,  Otsu())             

    HCa_binary_matrix = map(p -> p > 0.5 ? 1 : 0, HCa)
    On_binary_matrix = map(p -> p > 0.5 ? 1 : 0, On)

    HCa_On_binary_matrix = HCa_binary_matrix + On_binary_matrix

    out = fill(RGB{N0f8}(0,0,0), size(HCa))
    @inbounds for j in CartesianIndices(HCa_On_binary_matrix)
        if HCa_On_binary_matrix[j] > 1
            out[j] = RGB{N0f8}(0,1,0)     # green (overrides red beneath)
        elseif HCa_On_binary_matrix[j] == 1
            out[j] = RGB{N0f8}(1,0,0)     # red
        end
    end
    # saving generated image
    save(output_path, out)
end

function extract_sample_name(path::AbstractString)
    # Extract the filename from the path
    filename = split(path, '/')[end]
    # Remove the extension
    basename = split(filename, '.')[1]
    # Remove the channel (e.g., "-HCa" or "-On") before the number
    parts = split(basename, '-')
    # Assume the last part is the number, and the second to last is the channel
    # Remove the channel part
    sample_name = join(vcat(parts[1:end-2], parts[end]), "-")
    return sample_name
end

# Example usage
paths_HCa = readdir("./OsteonSurfaceAnalysis/DATA/FM40-1-R2/HCa"; join = true)
paths_On = readdir("./OsteonSurfaceAnalysis/DATA/FM40-1-R2/On"; join = true)
output_path = "./OsteonSurfaceAnalysis/DATA/FM40-1-R2/Processed_Images/test.png"
for (path_HCa, path_On) in zip(paths_HCa, paths_On)
    output_path = "./OsteonSurfaceAnalysis/DATA/FM40-1-R2/Processed_Images/" * extract_sample_name(path_HCa) * ".png"
    process_mask_data(path_HCa, path_On, output_path)
end


"""
make_tdelay_mask(Tdelay, path_to_processed_img, path_to_HCa, path_to_output; downsample=1)

Reads the processed (red/green) image and the HCa image, computes the Tdelay
zero-level contour, fills its interior, merges with the HCa mask, and writes a
color PNG:
- red  = HCa-only
- green = overlap of HCa and Tdelay interior

Arguments
---------
Tdelay :: Real
path_to_processed_img :: AbstractString   # e.g. ".../Processed_Images/...png"
path_to_HCa           :: AbstractString   # e.g. ".../HCa/...bmp"
path_to_output        :: AbstractString   # e.g. ".../mask0000.png"

Keyword
-------
downsample :: Int = 1
"""
function make_tdelay_mask(Tdelay::Real,
                          path_to_processed_img::AbstractString,
                          path_to_HCa::AbstractString,
                          path_to_output::AbstractString; downsample::Int=1)

    # --- constants ---
    BLACK = RGB{N0f8}(0, 0, 0)
    RED   = RGB{N0f8}(1, 0, 0)
    GREEN = RGB{N0f8}(0, 1, 0)

    # --- helpers ---
    edt(mask::BitArray) = sqrt.(DistanceTransforms.transform(boolean_indicator(mask)))
    edt_S(mask::BitArray) = edt(mask) .- edt(.!mask)

    function fill_internal_zeros(M::AbstractMatrix{<:Real})
        # Flood-fill zeros connected to the border, then flip the rest to 1
        filled = falses(size(M))
        queue = Tuple{Int,Int}[]
        H, W = size(M)

        # seed with all border zeros
        for i in 1:H, j in (1,W)
            if M[i,j] == 0 && !filled[i,j]
                push!(queue, (i,j)); filled[i,j] = true
            end
        end
        for j in 1:W, i in (1,H)
            if M[i,j] == 0 && !filled[i,j]
                push!(queue, (i,j)); filled[i,j] = true
            end
        end

        while !isempty(queue)
            i, j = popfirst!(queue)
            for (di,dj) in ((1,0),(-1,0),(0,1),(0,-1))
                ni, nj = i+di, j+dj
                if 1 ≤ ni ≤ H && 1 ≤ nj ≤ W && M[ni,nj] == 0 && !filled[ni,nj]
                    filled[ni,nj] = true
                    push!(queue, (ni,nj))
                end
            end
        end

        R = copy(M)
        for i in 1:H, j in 1:W
            if M[i,j] == 0 && !filled[i,j]
                R[i,j] = 1
            end
        end
        return R
    end

    # --- load processed (red/green) image and build masks ---
    img = load(path_to_processed_img)
    img = reverse(img, dims=2)
    img = img[1:downsample:end, 1:downsample:end]

    m_green = (img .== GREEN)
    m_red   = (img .== RED)

    H, W = size(img)
    Z = 1

    outer = falses(H, W, Z)    # region outside cement wall (red or green)
    inner = trues(H, W, Z)     # complement of green

    outer[:, :, 1] .= m_green .| m_red
    inner[:, :, 1] .= .!m_green

    # --- signed distance fields ---
    outer_dt_S = edt_S(outer)
    inner_dt_S = edt_S(inner)

    ϕ = (1 - Tdelay) .* outer_dt_S .- (Tdelay) .* inner_dt_S  # level-set field

    # --- zero contour → rasterize to matrix, then fill interior ---
    H, W = size(ϕ)[1:2]
    x = collect(1:H)
    y = collect(1:W)
    cset = CTR.contours(x,y,ϕ[:, :, 1], [0])
    if isempty(CTR.levels(cset)) || isempty(CTR.lines(first(CTR.levels(cset))))
        error("No zero-level contour found for the given Tdelay.")
    end
    line = first(CTR.lines(first(CTR.levels(cset))))
    xs, ys = CTR.coordinates(line)

    output_Matrix = zeros(Int, H, W)
    @inbounds for k in eachindex(xs)
        xi = round(Int, xs[k])
        yi = round(Int, ys[k])
        if 1 ≤ xi ≤ H && 1 ≤ yi ≤ W
            output_Matrix[xi, yi] = 1
        end
    end
    filled_Tdelay_matrix = fill_internal_zeros(output_Matrix)

    # --- build HCa mask (Otsu on blurred HCa image) ---
    hca = load(path_to_HCa)
    reverse!(hca, dims=2)
    hca_blur = imfilter(hca, Kernel.gaussian(2))
    HCa = binarize(hca_blur, Otsu())
    HCa_binary_matrix = map(p -> p > 0.5 ? 1 : 0, Gray.(HCa))

    # --- merge and colorize ---
    sum_mask = HCa_binary_matrix .+ filled_Tdelay_matrix

    out = fill(RGB{N0f8}(0,0,0), size(HCa_binary_matrix))
    @inbounds for j in CartesianIndices(sum_mask)
        if sum_mask[j] > 1
            out[j] = RGB{N0f8}(0,1,0)   # green = overlap
        elseif sum_mask[j] == 1
            out[j] = RGB{N0f8}(1,0,0)   # red = HCa only (or Tdelay only)
        end
    end

    save(path_to_output, out)
    return path_to_output
end

# Example usage
Tdelay = 0.1
path_to_processed_img = readdir("./OsteonSurfaceAnalysis/DATA/FM40-1-R2/Processed_Images"; join=true)[end] 
path_to_HCa = readdir("./OsteonSurfaceAnalysis/DATA/FM40-1-R2/HCa"; join = true)[end]
path_to_output = "./test2.png"

make_tdelay_mask(Tdelay,path_to_processed_img,path_to_HCa,path_to_output)

