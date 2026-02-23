using Plots
using Tables, CSV
include("network-generator.jl")
include("paths.jl")

"""
    count_flips(tracking::AbstractVector)

TBW
"""
function count_flips(tracking::AbstractVector)
    l = length(tracking)
    N = length(tracking[1])
    T = zeros(Int, N)

    for i ∈ 1:l-1, j ∈ 1:N
        T[j] += abs(tracking[i][j] - tracking[i+1][j])
    end

    count = 0
    # count0 = 0
    for flips ∈ T
        flips > 1 ? count += 1 : nothing
    end
    # state = tracking[1]
    # for i ∈ eachindex(T)
    #     if T[i] > 1
    #         count += 1
    #         if state[i] == 0
    #             count0 += 1
    #         end
    #     end
    # end

    return count / N# , count0 / count
end

"""
    count_ffrac(tracking::AbstractVector)

TBW
"""
function count_ffrac(tracking::AbstractVector)
    state = tracking[end]
    n = [0, 0]
    for value ∈ state
        n[value+1] += 1
    end
    return n[1] / length(state)
end

"""
    make_graphs(
    N::Int=N,
    p::Float64=p,
    num_realizations::Int=num_realizations
)

TBW
"""
function make_graphs(
    N::Int=N,
    p::Float64=p,
    num_realizations::Int=num_realizations
)

    graph_list = Graph[]

    for i ∈ 1:num_realizations
        g = erdos_renyi(N, p)
        @info "Made new ER graph, realization $(i)/$(num_realizations)"
        push!(graph_list, g)
    end

    return graph_list
end

"""
    results_plotter(
    sol_matrix::AbstractMatrix,
    x_vec::AbstractVector;
    ylabel::AbstractString,
    title::AbstractString,
    labels::AbstractVector
)

TBW
"""
function results_plotter(
    sol_matrix::AbstractMatrix,
    x_vec::AbstractVector;
    ylabel::AbstractString,
    title::AbstractString,
    labels::AbstractVector
)
    _, l_F = Base.size(sol_matrix)

    plt = plot(palette=:dracula,
        xlabel="θ", ylabel=ylabel,
        title=title,
        linewidth=0.75)

    for i ∈ 1:l_F
        plot!(plt, x_vec, sol_matrix[:, i],
            label=labels[i])
    end

    return plt
end

"""
    write_to_file(
    result::AbstractMatrix;
    directory::AbstractString,
    file_name::AbstractString,
    row_names::AbstractVector,
    col_names::AbstractVector
)

TBW
"""
function write_to_file(
    result::AbstractMatrix;
    directory::AbstractString,
    file_name::AbstractString,
    row_names::AbstractVector,
    col_names::AbstractVector
)

    file = directory *
           file_name * " " *
           join(string.(row_names), ", ") *
           ".csv"

    table = Tables.table(result)
    header = string.(col_names)

    CSV.write(file, table, header=header)
    @info "Wrote $(file_name) CSV"

    return nothing
end

"""
    giant_component_by_opinion(
    g::Graph,
    state::AbstractVector
)

Number of nodes in the GC of the induced sub-graph of nodes with `opinion`
"""
function giant_component_by_opinion(
    g::Graph,
    state::AbstractVector,
    op::Int,
)

    isg_edges = Tuple[]
    for (u, v) ∈ edges(g)
        if state[u] == op && state[v] == op
            push!(isg_edges, (u, v))
        end
    end

    if isempty(isg_edges)
        return 0
    else
        isg = graph_from_edges(isg_edges)
        gc = giant_component(isg)

        return length(gc)
    end
end
