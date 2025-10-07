import OsteocyteSpatDistribution as OSD

pcf_R = 95
pcf_ΔR = 0.5

# circle : Uniform
# generate osteocytes
N = 30; R = collect(range(10.0, 30.0, length=3)); distribution="Uniform"; direction="outward"
circle_Osteocytes = OSD.generate_osteocyte_points(N, R, distribution);

ρ = N / (π * R[end]^2)  # density of osteocytes
_circle_pcf_small, Radii= OSD.est_pcf2(circle_Osteocytes, π/6, pcf_R, pcf_ΔR; direction=direction)
_circle_pcf_inter, _= OSD.est_pcf2(circle_Osteocytes, π/3, pcf_R, pcf_ΔR; direction=direction)
_circle_pcf_large, _= OSD.est_pcf2(circle_Osteocytes, 3π/2, pcf_R, pcf_ΔR; direction=direction)


# Ellipse : Uniform
# generate osteocytes
distribution = "Uniform"; a = 1.5; b = 1.0;
ellipse_Osteocytes = OSD.generate_osteocyte_points(N, R, distribution; a, b);

ρ = N / (π * a * b * R[end]^2)  # density of osteocytes
_ellipse_pcf_small, _= OSD.est_pcf2(ellipse_Osteocytes, π/6, pcf_R, pcf_ΔR; direction=direction)
_ellipse_pcf_inter, _= OSD.est_pcf2(ellipse_Osteocytes, π/3, pcf_R, pcf_ΔR; direction=direction)
_ellipse_pcf_large, _= OSD.est_pcf2(ellipse_Osteocytes, 3π/2, pcf_R, pcf_ΔR; direction=direction)

# Circle : Random
# generate osteocytes
distribution = "Random";
rand_circle_Osteocytes = OSD.generate_osteocyte_points(N, R, distribution);

ρ = N / (π * R[end]^2)  # density of osteocytes
_rand_circle_pcf_small, _= OSD.est_pcf2(rand_circle_Osteocytes, π/6, pcf_R, pcf_ΔR; direction=direction)
_rand_circle_pcf_inter, _= OSD.est_pcf2(rand_circle_Osteocytes, π/3, pcf_R, pcf_ΔR; direction=direction)
_rand_circle_pcf_large, _= OSD.est_pcf2(rand_circle_Osteocytes, 3π/2, pcf_R, pcf_ΔR; direction=direction)



using CairoMakie

f_1 = OSD.plot_estpcf_for_idx(circle_Osteocytes, ellipse_Osteocytes, rand_circle_Osteocytes, 
                             _circle_pcf_small, _ellipse_pcf_small, _rand_circle_pcf_small,
                             _circle_pcf_inter, _ellipse_pcf_inter, _rand_circle_pcf_inter,
                             _circle_pcf_large, _ellipse_pcf_large, _rand_circle_pcf_large,
                             R, Radii, a, b, N, 1)

f_30 = OSD.plot_estpcf_for_idx(circle_Osteocytes, ellipse_Osteocytes, rand_circle_Osteocytes, 
                             _circle_pcf_small, _ellipse_pcf_small, _rand_circle_pcf_small,
                             _circle_pcf_inter, _ellipse_pcf_inter, _rand_circle_pcf_inter,
                             _circle_pcf_large, _ellipse_pcf_large, _rand_circle_pcf_large,
                             R, Radii, a, b, N, 30)

# Save the figures
save("Synthetic/Figures/more_points_pcf_est_outward_example_2.pdf", f_1)







"""
# All PCF figure
f = Figure(size = (1200, 1600), fontsize=22)
# --- Left plot: Circle points ---
ax_circle_top = Axis(
    f[1, 1],
    aspect=1,
    title = "Circle (Uniform)",
    xlabel = "x (μm)",
    ylabel = "y (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
ax_ellipse_top = Axis(
    f[1, 2],
    aspect=1,
    title = "Ellipse (Uniform)",
    xlabel = "x (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
ax_random_top = Axis(
    f[1, 3],
    aspect=1,
    title = "Circle (Random)",
    xlabel = "x (μm)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
# Plot concentric rings (radii)
for r in R
    θ = range(0, 2π, length=200)
    x = r * cos.(θ)
    y = r * sin.(θ)
    lines!(ax_circle_top, x, y, color=:gray, alpha=0.3, linewidth=1.5)
    lines!(ax_random_top, x, y, color=:gray, alpha=0.3, linewidth=1.5)
end 

for r in R
    θ = range(0, 2π, length=200)
    x = r * a * cos.(θ)
    y = r * b* sin.(θ)
    lines!(ax_ellipse_top, x, y, color=:gray, alpha=0.3, linewidth=1.5)
end 
# Plot random points on the circle
for (osteocyte_circ, osteocyte_ellipse, osteocyte_rand) in zip(circle_Osteocytes, ellipse_Osteocytes, rand_circle_Osteocytes)
    points_circ = osteocyte_circ.position
    scatter!(ax_circle_top, points_circ[1], points_circ[2], color=:dodgerblue, markersize=10)
    points_ellipse = osteocyte_ellipse.position
    scatter!(ax_ellipse_top, points_ellipse[1], points_ellipse[2], color=:dodgerblue, markersize=10)
    points_rand = osteocyte_rand.position
    scatter!(ax_random_top, points_rand[1], points_rand[2], color=:dodgerblue, markersize=10)
end

ax_circle_pcf_low = Axis(
    f[2, 1],
    aspect=1,
    ylabel = "g(r)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)

ax_ellipse_pcf_low = Axis(
    f[2, 2],
    aspect=1,
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)

ax_rand_pcf_low = Axis(
    f[2, 3],
    aspect=1,
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)

for (pcf_circ, pcf_ellipse, pcf_rand) in zip(_circle_pcf_small, _ellipse_pcf_small, _rand_circle_pcf_small)
    lines!(ax_circle_pcf_low, Radii, pcf_circ ./ N, linewidth=3);
    lines!(ax_ellipse_pcf_low, Radii, pcf_ellipse ./ N, linewidth=3);
    lines!(ax_rand_pcf_low, Radii, pcf_rand ./ N, linewidth=3);
end
ax_circle_pcf_inter = Axis(
    f[3, 1],
    aspect=1,
    ylabel = "g(r)",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
ax_ellipse_pcf_inter = Axis(
    f[3, 2],
    aspect=1,
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
ax_rand_pcf_inter = Axis(
    f[3, 3],
    aspect=1,
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
for (pcf_circ, pcf_ellipse, pcf_rand) in zip(_circle_pcf_inter, _ellipse_pcf_inter, _rand_circle_pcf_inter)
    lines!(ax_circle_pcf_inter, Radii, pcf_circ ./ N, linewidth=3);
    lines!(ax_ellipse_pcf_inter, Radii, pcf_ellipse ./ N, linewidth=3);
    lines!(ax_rand_pcf_inter, Radii, pcf_rand ./ N, linewidth=3);
end
ax_circle_pcf_large = Axis(
    f[4, 1],
    aspect=1,
    ylabel = "g(r)",
    xlabel = "r",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
ax_ellipse_pcf_large = Axis(
    f[4, 2],
    aspect=1,
    xlabel = "r",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
ax_rand_pcf_large = Axis(
    f[4, 3],
    aspect=1,
    xlabel = "r",
    xticklabelsize = 18,
    yticklabelsize = 18,
    xgridvisible = false,
    ygridvisible = false,
    backgroundcolor = :white
)
for (pcf_circ, pcf_ellipse, pcf_rand) in zip(_circle_pcf_large, _ellipse_pcf_large, _rand_circle_pcf_large)
    lines!(ax_circle_pcf_large, Radii, pcf_circ ./ N, linewidth=3);
    lines!(ax_ellipse_pcf_large, Radii, pcf_ellipse ./ N, linewidth=3);
    lines!(ax_rand_pcf_large, Radii, pcf_rand ./ N, linewidth=3);
end
f
"""

