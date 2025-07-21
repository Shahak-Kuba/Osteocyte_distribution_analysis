
# generate osteocytes
N = 10; R = range(10.0, 50.0, length=5); distribution = "Uniform";
Osteocytes = generate_osteocyte_points(N,R, distribution);

using CairoMakie

f = Figure(size = (400, 400), fontsize=22)

# --- Left plot: Circle points ---
ax = Axis(
    f[1, 1],
    aspect = 1,
    title = "Circle (Uniform)",
    xlabel = "x (μm)",
    ylabel = "y (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
# Plot concentric circles (radii)
for r in R
    θ = range(0, 2π, length=200)
    x = r * cos.(θ)
    y = r * sin.(θ)
    lines!(ax, x, y, color=:gray, alpha=0.3, linewidth=1.5)
end 
# Plot random points on the circle
for osteocyte in Osteocytes
    points = osteocyte.position
    scatter!(ax, points[1], points[2], color=:dodgerblue, markersize=10)
end
f

test_pcf= est_pcf2(Osteocytes, π/4, 100, 20)
