include("simulate.jl")
using DataStructures

function radius(g::Graph, k::Vector{Int}, θ::Float64, node::Int64, state::Vector{Int})
    N = size(g)
    visited = falses(N)
    dists::Vector{Int} = fill(-1, N)
    affected = Set{Int64}()

    Q = Queue{Int64}()
    enqueue!(Q, node)
    visited[node] = true
    dists[node] = 0

    current_state = Base.copy(state)

    current_state[node] = Int(1 - current_state[node])

    while !isempty(Q)
        current = dequeue!(Q)

        for nbr ∈ neighbors(g, current)
            if !visited[nbr]
                new = new_state(g, k, θ, node, current_state)

                if new ≠ current_state[nbr]
                    current_state[nbr] = new
                    visited[nbr] = true
                    dists[nbr] = dists[current] + 1
                    enqueue!(Q, nbr)
                    push!(affected, nbr)
                end
            end
        end
    end

    return maximum(dists), affected
end

function calculate_instability(g::Graph, θ::Float64, state::Vector{Int})
    N = size(g)
    k = degrees(g)
    
    radii = zeros(Int, N)
    instability = zeros(Int, N) # node wise stability (larger is more unstable)
    for node ∈ nodes(g)
        r, aff = radius(g, k, θ, node, state)
        radii[node] = r
        for affected ∈ aff
            instability[affected] += 1
        end
    end

    c = countmap(instability) # counts the number of nodes with a particular stability value

    p_stab = Dict{Int,Float64}()
    for key ∈ keys(c)
        p_stab[key] = c[key] / N
    end

    return p_stab
end

function instability_mode(p_stab::AbstractDict)
    max_key = reduce((x, y) -> p_stab[x] ≥ p_stab[y] ? x : y, keys(p_stab))

    return p_stab[max_key]
end

function instability_avg(p_stab::AbstractDict)
    sum = 0

    for key in keys(p_stab)
	    sum += p_stab[key]
    end

    return sum / length(keys(p_stab))
end
