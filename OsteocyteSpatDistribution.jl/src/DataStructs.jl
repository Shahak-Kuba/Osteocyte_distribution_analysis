struct Osteocyte
    position::Tuple{Float64, Float64}
    dist_from_origin::Float64
    τ::Tuple{Float64, Float64}  # Orientation vector
    n::Tuple{Float64, Float64} # Normal vector
end