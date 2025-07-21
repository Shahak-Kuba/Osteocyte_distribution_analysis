struct Ellipse
    center::Vector{Float64}
    axes::Vector{Float64}  # [a, b] semi-axes lengths
    angle::Float64         # orientation angle in radians
end