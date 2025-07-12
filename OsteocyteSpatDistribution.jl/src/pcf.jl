function GaussianKernel(x)
    return (1 ./ (2 .* π) ) .* exp.(-0.5 .* x.^2)
end

function λ_est(points, Kernel, h)
    """
    Estimate the intensity function λ at each point using a kernel density estimator.

    Args:
        points: A vector of points (x, y) where the intensity is estimated.
        Kernel: The kernel function to use for estimation.
        h: The bandwidth parameter for the kernel.

    Returns:
        A vector of estimated intensities at each point.
    """
    n = length(points)
    λ = zeros(n)

    for i in 1:n
        xi, yi = points[i]
        sum_kernel = 0.0

        for j in 1:n
            xj, yj = points[j]
            distance = sqrt((xi - xj)^2 + (yi - yj)^2)
            sum_kernel += Kernel(distance / h)
        end

        λ[i] = sum_kernel / (n * h^2)
    end

    return λ
end