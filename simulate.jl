include("lib/graph.jl")
using StatsBase, Random

function new_state(
    g::Graph,
    k::Vector{Int},
    θ::Float64,
    node::Int64,
    state::Vector{Int})

    m0::Int = 0
    m1::Int = 0

    for nbr ∈ neighbors(g, node)
        state[nbr] == 1 ? m1 += k[nbr] : m0 += k[nbr]
    end

    m0 + m1 > 0 ? node_frac = m1 / (m0 + m1) : node_frac = 0

    node_frac ≥ θ ? new_state = 1 : node_frac < θ ? new_state = 0 : new_state = state[node]

    return new_state
end

function fractional_IMR(
    g::Graph;
    θ::Float64,
    F₀::Float64,
    max_sweeps::Int=100)

    k = degrees(g)
    N = size(g)
    state::Vector{Int} = sample([0, 1], Weights([F₀, 1 - F₀]), N, replace=true)

    full_tracking = Vector{Int64}[]
    full_tracking = push!(full_tracking, state)

    for _ ∈ 1:max_sweeps
        changes = 0
        current = Base.copy(full_tracking[end])

        for node ∈ shuffle(collect(nodes(g)))
            new = new_state(g, k, θ, node, current)

            if new ≠ current[node]
                current[node] = new
                changes += 1
            end
        end
        
        full_tracking = push!(full_tracking, current)
        
        if changes == 0
            break
        end
    end

    return full_tracking
end
