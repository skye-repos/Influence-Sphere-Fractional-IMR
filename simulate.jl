include("lib/graph.jl")
using StatsBase, Random

"""
    new_state(g::Graph, k::AbstractVector, θ::Float64, node::Integer, state::AbstractVector)

Compute the new state of `node` in `g` from the given `state` based on some logic using `θ`.
"""
function new_state(
    g::Graph,
    k::AbstractVector,
    θ::Float64,
    node::Integer,
    state::AbstractVector
)
    # Set-up vars for node fraction calculations
    m0::Int = 0
    m1::Int = 0

    # Actual Logic for computing the node fraction
    for nbr ∈ neighbors(g, node)
        state[nbr] == 1 ? m1 += k[nbr] : m0 += k[nbr]
    end

    m0 + m1 > 0 ? frac = m1 / (m0 + m1) : frac = 0

    # Computing the new_state based on an algorithm
    frac ≥ θ ? new_state = 1 : frac < θ ? new_state = 0 : new_state = state[node]

    return new_state
end

"""
    fractional_IMR(g::Graph; θ::Float64, F₀::Float64, max_sweeps::Int=100)

Seed a fraction of nodes `F₀` in network `g` with opinion 0 and the rest with opinion 1, then loop over each node in a random order to compute the new opinions. Iterate until the network's opinions do not change at all.
"""
function fractional_IMR(
    g::Graph;
    θ::Float64,
    F₀::Float64,
    max_sweeps::Int=100
)

    k = degrees(g)
    # Seed F₀ * N nodes in `g` with 0 states
    state::Vector{Int} = sample([0, 1],
        Weights([F₀, 1 - F₀]),
        nv(g), replace=true)

    tracking = Vector{Int64}[]
    push!(tracking, state)

    for _ ∈ 1:max_sweeps
        changes = 0
        current = Base.copy(tracking[end])

        for node ∈ shuffle(collect(nodes(g)))
            new = new_state(g, k, θ, node, current)

            if new ≠ current[node]
                current[node] = new
                changes += 1
            end
        end

        push!(tracking, current)

        if changes == 0
            break
        end
    end

    return tracking
end
