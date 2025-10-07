function plot_estpcf_for_idx(circle_Osteocytes, ellipse_Osteocytes, rand_circle_Osteocytes, 
                             _circle_pcf_small, _ellipse_pcf_small, _rand_circle_pcf_small,
                             _circle_pcf_inter, _ellipse_pcf_inter, _rand_circle_pcf_inter,
                             _circle_pcf_large, _ellipse_pcf_large, _rand_circle_pcf_large,
                             R, Radii, a, b, N, idx)
    # Single PCF figure
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
        limits!(ax_ellipse_top, -maximum([a,b])*R[end]*1.1, maximum([a,b])*R[end]*1.1, -maximum([a,b])*R[end]*1.1, maximum([a,b])*R[end]*1.1)
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
    scatter!(ax_circle_top, circle_Osteocytes[idx].position[1], circle_Osteocytes[idx].position[2], color=:crimson, markersize=15)
    scatter!(ax_ellipse_top, ellipse_Osteocytes[idx].position[1], ellipse_Osteocytes[idx].position[2], color=:crimson, markersize=15)
    scatter!(ax_random_top, rand_circle_Osteocytes[idx].position[1], rand_circle_Osteocytes[idx].position[2], color=:crimson, markersize=15)

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
    lines!(ax_circle_pcf_low, Radii, _circle_pcf_small[idx] ./ N, linewidth=3, color=:crimson);
    lines!(ax_ellipse_pcf_low, Radii, _ellipse_pcf_small[idx] ./ N, linewidth=3, color=:crimson);
    lines!(ax_rand_pcf_low, Radii, _rand_circle_pcf_small[idx] ./ N, linewidth=3, color=:crimson);

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
    lines!(ax_circle_pcf_inter, Radii, _circle_pcf_inter[idx] ./ N, linewidth=3, color=:crimson);
    lines!(ax_ellipse_pcf_inter, Radii, _ellipse_pcf_inter[idx] ./ N, linewidth=3, color=:crimson);
    lines!(ax_rand_pcf_inter, Radii, _rand_circle_pcf_inter[idx] ./ N, linewidth=3, color=:crimson);

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
    lines!(ax_circle_pcf_large, Radii, _circle_pcf_large[idx] ./ N, linewidth=3, color=:crimson);
    lines!(ax_ellipse_pcf_large, Radii, _ellipse_pcf_large[idx] ./ N, linewidth=3, color=:crimson);
    lines!(ax_rand_pcf_large, Radii, _rand_circle_pcf_large[idx] ./ N, linewidth=3, color=:crimson);

    return f
end