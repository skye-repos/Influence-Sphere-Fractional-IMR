include("graph.jl")

"""
Computes the edge list for a connected graph with N0 nodes
"""
function connected_graph(N0::Int)
    adjacency_list::Vector{Vector{Int}} = []
    nodes = collect(1:N0)
    for i = 1:N0
        push!(adjacency_list, filter(e->e≠i, nodes))
    end
    return Graph(adjacency_list)
end

"""
Computes the next barabasi-albert iteration of a network `g` where a new node brings with it `m` half links.
"""
function ba_step!(g::Graph, m::Int64)
    # Create a pool of nodes weighted by their degree
    pool = vcat([fill(i, degrees(g)[i]) for i in 1:size(g)]...)
    targets = Int[]
    while length(targets) < m
        v = rand(pool)
        v ∉ targets && push!(targets, v)
    end

    for v in targets
        push!(g.adjacency_list[v], size(g)+1)
    end

    push!(g.adjacency_list, targets)

    return g
end

"""
Given an initial network `g`, computes the Barabasi-Albert graph after `n_steps` nodes are added each of which brings with it `m` half-links
"""
function barabasi_albert!(g::Graph, m::Int64, n_steps::Int64)
    for _ = 1:n_steps
        ba_step!(g, m)
    end
    return g
end

"""
Generate Erdos-Renyi random network G(N, p)
""" 
function erdos_renyi(N, p)
    adjacency_list = [Vector{Int64}() for _ = 1:N]
    g = Graph(adjacency_list)

    n_edges = Int(ceil(N*(N-1)*p/2))

    if p < 0.1
        for _ = 1:n_edges
            i = Int.(rand(1:N))
            j = Int.(rand(1:N))
            while j ∈ adjacency_list[i]
                j = Int.(rand(1:N))
            end
            add_link!(g, Int(i), Int(j))
        end
    else
        for i = 1:N
            for j = 1:N
                u = rand()
                if u < p
                    add_link!(g, Int(i), Int(j))
                end
            end
        end
    end
    return g
end

"""
Generate a ζ-model network of side length `L`, characteristic link-length of
`ζ`, and average degree `k`
"""
function ζ_model(L::Integer, ζ::Number, k_avg::Number)
    local N = L^2

    coords = Vector{Tuple{Int, Int}}()

    for i ∈ 1:L, j ∈ 1:L
        push!(coords, (i, j))
    end

    @inline function distance(a::Tuple{Int, Int}, b::Tuple{Int, Int})
        x = min(abs(a[1] - b[1]), L - abs(a[1] - b[1]))
        y = min(abs(a[2] - b[2]), L - abs(a[2] - b[2]))
        return sqrt(x^2 + y^2)
    end

    P = zeros(N, N)
    @inbounds for j ∈ 1:N, i ∈ 1:j-1
        P[i, j] = exp(-distance(coords[i], coords[j]) / ζ)
    end

    P .*= (N * k_avg / 2) / sum(P)

    g = Graph([Vector{Int64}() for _ ∈ 1:N])

    @inbounds for j ∈ 1:N, i ∈ 1:j-1
        rand() < P[i, j] && add_link!(g, Int(i), Int(j))
    end

    return g
end
