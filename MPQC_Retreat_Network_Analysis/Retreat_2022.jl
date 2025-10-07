# ========== Setup ==========
# If needed, add packages once:
# using Pkg; Pkg.add(["Graphs", "SimpleWeightedGraphs", "Statistics"])

using Graphs
using SimpleWeightedGraphs
using Statistics
using GraphMakie
using NetworkLayout

"""
    build_graph_and_stats(A::AbstractMatrix{<:Real}; tol=1e-8)

Given an interaction matrix `A` with values:
- 0.1 => no interaction (no edge)
- 0.5 => interaction (weight = 1.0)
- 1.0 => strong interaction (weight = 2.0)
Diagonal entries are ignored (assumed to be 1 for self-interaction).

Returns:
- `g`: a weighted, undirected graph
- `stats`: a NamedTuple with network statistics
"""
function build_graph_and_stats(A::AbstractMatrix{<:Real}; tol=1e-8)
    @assert size(A,1) == size(A,2) "Matrix A must be square"
    n = size(A,1)

    # Make sure the matrix is symmetric
    A_sym = max.(A, A')

    g = SimpleWeightedGraph(n)

    for i in 1:n-1, j in i+1:n
        # Skip diagonal (just to be safe, though we don't iterate over it)
        a = A_sym[i,j]
        if isapprox(a, 0.5; atol=tol)
            add_edge!(g, i, j, 1.0)
        elseif isapprox(a, 1.0; atol=tol)
            add_edge!(g, i, j, 2.0)
        else
            # If it's 0.1 or anything close → no edge
        end
    end

    # ---- Metrics ----
    deg = degree(g)
    local_cc = [local_clustering_coefficient(g, v) for v in 1:nv(g)]
    global_cc = global_clustering_coefficient(g)
    dens = Graphs.density(g)
    avg_deg = mean(deg)

    stats = (
        degrees = deg,
        local_clustering = local_cc,
        global_clustering = global_cc,
        density = dens,
        average_degree = avg_deg,
        nv = nv(g),
        ne = ne(g)
    )

    return g, stats
end

# using Pkg; Pkg.add(["Graphs", "SimpleWeightedGraphs", "Statistics"])

using Graphs
using SimpleWeightedGraphs
using Statistics

# using Pkg; Pkg.add(["Graphs", "SimpleWeightedGraphs", "Statistics"])
using Graphs
using SimpleWeightedGraphs
using Statistics

"""
    build_graph_and_stats(A::AbstractMatrix{<:Real}; tol=1e-8)

Build an undirected, weighted graph from an interaction matrix `A` with values:
- 0.1 => no edge
- 0.5 => edge (weight = 1.0)
- 1.0 => strong edge (weight = 2.0)
Diagonal is ignored (assumed self-interaction = 1.0).

Returns `(g, stats)` where `g` is a `SimpleWeightedGraph` and `stats` is a NamedTuple:

Basic:
- degrees::Vector{Int}
- local_clustering::Vector{Float64}
- global_clustering::Float64
- density::Float64
- average_degree::Float64
- nv::Int, ne::Int

Connectivity / growth:
- n_components::Int
- largest_component_size::Int
- average_path_length_lcc::Union{Float64,Missing}   # unweighted, on LCC
- diameter_lcc::Union{Float64,Missing}              # unweighted, on LCC
- degree_assortativity::Union{Float64,Missing}      # Newman r (unweighted)
"""
function build_graph_and_stats2(A::AbstractMatrix{<:Real}; tol=1e-8)
    @assert size(A,1) == size(A,2) "Matrix A must be square"
    n = size(A,1)

    # Defensively symmetrize (keep stronger of A[i,j], A[j,i]); diagonal ignored
    A_sym = max.(A, A')

    # Build weighted, undirected graph
    g = SimpleWeightedGraph(n)
    for i in 1:n-1, j in i+1:n
        a = A_sym[i,j]
        if isapprox(a, 0.5; atol=tol)
            add_edge!(g, i, j, 1.0)
        elseif isapprox(a, 1.0; atol=tol)
            add_edge!(g, i, j, 2.0)
        end
        # 0.1 (or other) => no edge
    end

    # --- Basic metrics (Graphs.jl topology is unweighted by default) ---
    deg = degree(g)
    local_cc = [local_clustering_coefficient(g, v) for v in 1:nv(g)]
    global_cc = global_clustering_coefficient(g)
    dens = Graphs.density(g)                 # qualify to avoid clashes
    avg_deg = mean(deg)

    # --- Connected components (unweighted semantics) ---
    comps = connected_components(g)          # Vector{Vector{Int}}
    ncomp = length(comps)
    lcc_size = maximum(length.(comps); init=0)

    # --- LCC path metrics (average path length & diameter) ---
    apl_lcc = missing
    diam_lcc = missing
    if lcc_size > 1
        # pick largest component
        lcc_idx = argmax(length.(comps))
        lcc = comps[lcc_idx]
        # map original vertex ids -> 1:|LCC|
        idxmap = Dict(v => i for (i, v) in enumerate(lcc))

        # build an unweighted SimpleGraph induced by LCC
        ug_lcc = SimpleGraph(length(lcc))
        for e in edges(g)
            u, v = src(e), dst(e)
            if haskey(idxmap, u) && haskey(idxmap, v)
                add_edge!(ug_lcc, idxmap[u], idxmap[v])
            end
        end

        # all-pairs shortest paths (unweighted => each edge length 1)
        fw = floyd_warshall_shortest_paths(ug_lcc)
        D = fw.dists

        finite_d = Float64[]
        diam = 0.0
        for i in 1:size(D,1), j in 1:size(D,2)
            if i != j && isfinite(D[i,j])
                push!(finite_d, D[i,j])
                if D[i,j] > diam
                    diam = D[i,j]
                end
            end
        end
        apl_lcc = isempty(finite_d) ? missing : mean(finite_d)
        diam_lcc = isempty(finite_d) ? missing : diam
    end

    # --- Degree assortativity (Newman r) on unweighted topology ---
    assort = missing
    if ne(g) > 0
        m = ne(g)
        s1 = 0.0   # sum(j*k) over edges
        s2 = 0.0   # sum((j+k)/2)
        s3 = 0.0   # sum((j^2 + k^2)/2)
        for e in edges(g)
            u, v = src(e), dst(e)
            ju, jv = deg[u], deg[v]
            s1 += ju * jv
            s2 += (ju + jv) / 2
            s3 += (ju^2 + jv^2) / 2
        end
        Ex_jk = s1 / m
        Ex_q  = s2 / m
        Ex_q2 = s3 / m
        denom = (Ex_q2 - Ex_q^2)
        assort = iszero(denom) ? missing : (Ex_jk - Ex_q^2) / denom
    end

    stats = (
        degrees = deg,
        local_clustering = local_cc,
        global_clustering = global_cc,
        density = dens,
        average_degree = avg_deg,
        n_components = ncomp,
        largest_component_size = lcc_size,
        average_path_length_lcc = apl_lcc,
        diameter_lcc = diam_lcc,
        degree_assortativity = assort,
        nv = nv(g),
        ne = ne(g)
    )

    return g, stats
end

 
function plot_network(g::AbstractGraph; node_labels::Bool=false, savepath::Union{Nothing,String}=nothing)
    n = nv(g)
    deg = degree(g)

    # Node size based on degree
    dmin, dmax = minimum(deg), maximum(deg)
    nodesize = dmax == dmin ? fill(14.0, n) : 10 .+ 18 .* ((deg .- dmin) ./ (dmax - dmin))

    # Edge styling based on weight
    edge_color = Symbol[]
    edge_width = Float64[]
    for e in edges(g)
        u, v = src(e), dst(e)
        w = try
            Graphs.weight(g, u, v)
        catch
            1.0
        end
        if w > 1.0
            push!(edge_color, :black)
            push!(edge_width, 3.5)
        else
            push!(edge_color, :gray)
            push!(edge_width, 1.8)
        end
    end

    # ✅ Use spring_layout from NetworkLayout.jl
    layout = Stress() 

    CairoMakie.activate!()
    fig = Figure(resolution=(900, 700))
    ax  = CairoMakie.Axis(fig[1, 1], title="Interaction Network",
               xticksvisible=false, yticksvisible=false,
               xgridvisible=false, ygridvisible=false)

    graphplot!(ax, g;
        layout,
        node_size = nodesize,
        node_color = :dodgerblue,
        edge_color = edge_color,
        edge_width = edge_width,
        nlabels   = node_labels ? (1:n) : nothing,
        nlabels_align = (:center, :center),
        nlabels_color = :black,
        nlabels_textsize = 10
    )

    hidespines!(ax)
    if savepath !== nothing
        save(savepath, fig)
        @info "Saved figure to $savepath"
    end

    return fig
end

MPQC_Network_Theme_2022 = [
    1.0   1.0   0.1   1.0   0.1   0.1   0.1;
    1.0   1.0   0.5   0.1   0.1   0.1   0.1;
    0.1   1.0   1.0   0.5   0.5   0.1   0.1;
    1.0   0.5   1.0   1.0   1.0   1.0   0.5;
    0.1   0.1   0.5   1.0   1.0   0.1   1.0;
    0.1   0.1   0.1   1.0   0.1   1.0   0.1;
    0.1   0.1   0.5   0.5   1.0   0.1   1.0
]

MPQC_Network_Theme_2025 = [
    1.0   1.0   0.1   1.0   0.1   0.1   0.1   1.0   0.1   0.1   0.1   0.1   0.1   0.1   0.1;
    1.0   1.0   0.5   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.5   0.1   0.1   0.1   0.1;
    0.1   1.0   1.0   0.5   0.5   0.1   0.1   0.5   1.0   0.1   0.5   0.1   1.0   0.1   0.1;
    1.0   0.5   1.0   1.0   1.0   0.5   1.0   1.0   0.5   1.0   0.5   0.5   0.5   0.5   0.5;
    0.1   0.1   0.5   1.0   1.0   1.0   0.1   0.5   0.1   0.1   0.1   1.0   0.1   1.0   1.0;
    0.1   0.1   0.1   0.5   1.0   1.0   0.1   0.1   0.1   0.1   0.1   0.5   0.1   1.0   1.0;
    0.1   0.1   0.1   1.0   0.1   0.1   1.0   0.5   0.1   0.1   0.1   0.1   0.1   0.1   0.1;
    1.0   0.1   0.5   1.0   0.5   0.1   0.5   1.0   0.5   1.0   0.1   0.5   0.5   0.5   0.1;
    0.1   0.1   1.0   0.5   0.5   0.1   0.1   0.5   1.0   0.1   0.1   0.1   1.0   0.1   0.1;
    0.1   0.1   0.1   1.0   0.1   0.1   0.1   1.0   0.1   1.0   0.1   0.1   0.1   0.1   0.1;
    0.1   0.5   0.5   0.1   0.1   0.1   0.1   0.1   0.5   0.1   1.0   0.1   0.1   0.1   0.1;
    0.1   0.1   0.1   0.5   1.0   0.5   0.1   0.5   0.1   0.1   0.1   1.0   0.1   0.5   0.5;
    0.1   0.1   1.0   0.5   0.1   0.1   0.1   0.5   1.0   0.1   0.1   0.1   1.0   0.1   0.1;
    0.1   0.1   0.1   1.0   1.0   1.0   0.1   0.5   0.1   0.1   0.1   1.0   0.1   1.0   0.1;
    0.1   0.1   0.5   0.5   1.0   1.0   0.1   0.5   0.1   0.1   0.1   1.0   0.1   1.0   1.0
]

g_2022, stats_2022 = build_graph_and_stats2(MPQC_Network_Theme_2022)

g_2025, stats_2025 = build_graph_and_stats2(MPQC_Network_Theme_2025)

println("Nodes 2022: ", stats_2022.nv, " | Edges 2022: ", stats_2022.ne)
println("Nodes 2025: ", stats_2025.nv, " | Edges 2025: ", stats_2025.ne)

println("Degrees 2022: ", stats_2022.degrees)
println("Degrees 2025: ", stats_2025.degrees)

println("Local clustering coefficients 2022: ", stats_2022.local_clustering)
println("Local clustering coefficients 2025: ", stats_2025.local_clustering)

println("Global clustering coefficient 2022: ", stats_2022.global_clustering)
println("Global clustering coefficient 2025: ", stats_2025.global_clustering)

println("Density 2022: ", stats_2022.density)
println("Density 2025: ", stats_2025.density)

println("Average degree 2022: ", stats_2022.average_degree)
println("Average degree 2025: ", stats_2025.average_degree)

println("Average degree 2022: ", stats_2022.average_path_length_lcc)
println("Average degree 2025: ", stats_2025.average_path_length_lcc)


plot_network(g_2022)
plot_network(g_2025)
