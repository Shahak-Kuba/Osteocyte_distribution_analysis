import OsteocyteSpatDistribution as OSD
using CairoMakie
using Random

Random.seed!(45)  # For reproducibility
Distribution_type = "Uniform" 
N = 40
Radii = range(10.0, 50.0, length=5)
points = generate_osteocyte_points(Distribution_type, N, Radii)

f = Figure(size = (800, 400), fontsize=22)
ax = Axis(
    f[1, 1],
    aspect = 1,
    title = "Synthetic Osteocytes",
    xlabel = "x (μm)",
    ylabel = "y (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)   
# Plot concentric circles (radii)
for r in Radii
    θ = range(0, 2π, length=200)
    x = r * cos.(θ)
    y = r * sin.(θ)
    lines!(ax, x, y, color=:gray, alpha=0.3, linewidth=1.5)
end
# Plot osteocyte points
for (x, y) in points
    scatter!(ax, x, y, color=:dodgerblue, markersize=10, label="Osteocyte")
end
f