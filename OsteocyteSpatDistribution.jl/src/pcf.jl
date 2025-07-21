function est_pcf2(cell::Osteocyte, all_cells:: Vector{Osteocyte}, α, Rmax, ΔR; direction="inward")
    # cell: Osteocyte data struct containing cell information
    # α: angle in radians
    # Rmax: maximum radius for the PCF calculation
    # ΔR: radial bin width
    # direction: "inward" or "outward"
    
    # check that Rmax is divisible by ΔR
    @assert Rmax % ΔR == 0 "Rmax must be divisible by ΔR"

    # defining normalised directional vector
    if direction == "inward" 
        normal_dir = cell.n ./ norm(cell.n)
    elseif direction == "outward"
        normal_dir = -1 .* (cell.n ./ norm(cell.n))
    else
        throw(ArgumentError("direction must be 'inward' or 'outward'"))
    end
    Radii = range(0, stop=Rmax, step=ΔR)
    # Initialize the PCF array
    n_bins = Int(Rmax / ΔR)
    pcf = zeros(n_bins)

    # calculate pcf for each radius
    for other_cell in all_cells
        if other_cell != cell
            # Calculate the vector from the cell to the other cell
            vec_to_other = other_cell.position .- cell.position
            # Calculate the distance to the other cell
            dist = norm(vec_to_other)
            # Calculate angle between normal vector and vector to other cell
            angle = acos(clamp(dot(normal_dir, vec_to_other ./ dist), -1.0, 1.0))
            if angle < α/2
                bin_index = Int(floor(dist / ΔR)) + 1  # +1 for 1-based indexing in Julia
                if bin_index <= n_bins && bin_index > 0  # Ensure the index is within bounds
                    pcf[bin_index] += 1  # Increment the count for this bin
                end 
            end
        end
    end
    return pcf
end

function est_pcf2(cells::Vector{Osteocyte}, α, Rmax, ΔR; direction="inward")
    # cells: Vector of Osteocyte data structs
    # α: angle in radians
    # Rmax: maximum radius for the PCF calculation
    # ΔR: radial bin width
    # direction: "inward" or "outward"
    
    pcf_results = []
    for cell in cells
        pcf = est_pcf2(cell, cells, α, Rmax, ΔR; direction=direction)
        push!(pcf_results, pcf)
    end
    return pcf_results
end 

