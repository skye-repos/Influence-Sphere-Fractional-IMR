include("simulate.jl")
using DataStructures

"""
    influence_sphere(
    g::Graph,
    k::AbstractVector,
    θ::Float64,
    node::Int64,
    state::Vector{Int})

Creates a Directed Acyclic Graph of nodes that are within the 'pure' influence
sphere of `node` given the topology of `g` and initial `state`.
"""
function influence_sphere(
    g::Graph,
    k::AbstractVector,
    θ::Float64,
    node::Int64,
    state::Vector{Int}
)
    # Allocate some objects for the BFS
    visited = falses(nv(g))
    dists::Vector{Int} = fill(-1, nv(g))
    affected = Set{Int}()

    # Start the BFS Queue
    Q = Queue{Int}()
    enqueue!(Q, node)
    visited[node] = true
    dists[node] = 0

    current_state::Vector{Int} = Base.copy(state)
    flipped_state::Vector{Int} = Base.copy(state)
    flipped_state[node] = 1 - flipped_state[node]

    while !isempty(Q)
        current = dequeue!(Q)

        for nbr ∈ neighbors(g, current)
            if !visited[nbr]
                new_curr = new_state(g, k, θ, nbr,
                    current_state)
                new_flip = new_state(g, k, θ, nbr,
                    flipped_state)

                # Check if new state of a neighbor due to the flip is 1)
                # different from its original because of the flip 2) wouldn't
                # have changed even if the flip hadn't happened
                if new_flip ≠ current_state[nbr] && new_curr == current_state[nbr]
                    flipped_state[nbr] = new_flip
                    visited[nbr] = true
                    dists[nbr] = dists[current] + 1
                    enqueue!(Q, nbr)
                    push!(affected, nbr)
                end
            end
        end
    end
    return dists, affected
end

"""
    total_instability(g::Graph, θ::Float64, state::Vector{Int})

Compute the node-wise instability of `g` under the fractional influence majority
opinion dynamics using the influence-sphere for `θ`, with a given `state`.
"""
function total_instability(g::Graph, θ::Float64, state::Vector{Int})
    # Node wise instability is defined by the number of times a particular node
    # is affected/in the influence sphere of another node.  one definition is a
    # simple tally another is to modulate the contribution based on their
    # distance/path length of influence (current)
    N = nv(g)
    k = degrees(g)
    I = zeros(Float64, N)

    for node ∈ nodes(g) # loop over all nodes
        # Compute the IS & the Distances
        D, A = influence_sphere(g, k, θ, node, state)
        # For all a ∈ IS(node), add contribution to instability
        for a ∈ A
            I[a] += 1 / D[a]
        end
    end

    I = I / (N - 1)

    return I
end

"""
    opinion_instability(I::AbstractVector, state::AbstractVector)

Given the node-wise instability `I` corresponding to a given `state`, separate
out the instabilities of opinion 0 and opinion 1 nodes and compute their
expectation values.
"""
function opinion_instability(I::AbstractVector, state::AbstractVector)
    N = length(I)
    I0 = zeros(Float64, N)
    I1 = zeros(Float64, N)

    is0 = state .== 0
    N0 = count(is0)
    F0 = N0 / N
    N1 = N - N0
    F1 = 1 - F0
    for node ∈ eachindex(I)
        state[node] == 0 ? I0[node] = I[node] : I1[node] = I[node]
    end

    I0 = I0 / F0
    I1 = I1 / F1

    return I0, I1
end

"""
    expected_instability(V::AbstractVector, bin::Float64=1e-5)

TBW
"""
# function expected_instability(V::AbstractVector, bin::Float64)
#     if isempty(V)
#         return 0.0
#     end

#     Vₘ = maximum(V)
#     if Vₘ == 0
#         return 0.0
#     end

#     N = length(V)
#     v = range(0, Vₘ, step=bin)
#     p = zeros(Float64, length(v))

#     for ins ∈ V
#         i = floor(Int, ins / bin) + 1
#         p[i] += 1 / N
#     end

#     return v' * p
# end

function expected_instability(V::AbstractVector)
    if isempty(V)
        return 0.0
    end

    Vₘ = maximum(V)
    N = length(V)
    if Vₘ == 0
        return 0.0
    end

    return sum(V) / N
end
