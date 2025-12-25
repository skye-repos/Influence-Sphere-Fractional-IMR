using Plots
using Tables, CSV
include("lib/network-generator.jl")

function count_flips_fraction(full_tracking::Vector{Vector{Int}})
    l = length(full_tracking)
    N = length(full_tracking[1])
    tracker = zeros(Int, N)

    for i ∈ 1:l-1, j ∈ 1:N
        tracker[j] += abs(full_tracking[i][j] - full_tracking[i+1][j])
    end

    count = 0

    for flips ∈ tracker
        if flips > 1
            count += 1
        end
    end

    return count / N
end

function final_fraction(full_tracking)
    state = full_tracking[end]
    n = [0, 0]
    for value in state
        n[value+1] += 1
    end
    return n[1] / length(state)
end

function make_graphs(
    dir;
    type="ER",
    N=N, p=p,
    N0=N0, m=m,
    num_realizations=num_realizations)

    dir = dir * "/$(type)/N = $(N)/"

    graph_list = Vector{Graph}(undef, num_realizations)

    for i ∈ 1:num_realizations
        if type == "ER"
            name = dir * "p = $(p), realization #$(i).bin"
        elseif type == "BA"
            name = dir * "N0 = $(N0), m = $(m), realization #$(i).bin"
        else
            error("Unsupported graph kind not ER or BA")
        end

        if type == "ER"
            g = erdos_renyi(N, p)
            @info "Made new ER graph, realization $(i)/$(num_realizations)"
        elseif type == "BA"
            g = connected_graph(N0)
            barabasi_albert!(g, m, N - N0)
            @info "Made new BA graph, realization $(i)/$(num_realizations)"
        end

        graph_list[i] = g
    end

    return graph_list
end

function results_plotter(sol_matrix::AbstractMatrix, θ_list=θ_list; ylabel, title, labels=plot_labels)
    l_θ, l_F = Base.size(sol_matrix)

    plt = plot(palette=:dracula,
        xlabel="θ", ylabel=ylabel,
        title=title,
        linewidth=0.75)

    for i ∈ 1:l_F
        plot!(plt, θ_list, sol_matrix[:, i],
            label=label = labels[i])
    end

    return plt
end

"""
Convert a vector of vectors to a matrix
"""
function vecvec_to_matrix(results::Vector{Vector{Float64}})
    sol = vcat([], results[1])

    for n ∈ eachindex(results)[2:end]
        sol = hcat(sol, results[n])
    end

    return sol
end

"""
Write a matrix `result` to file with name `file_name` and header `header`
"""
function write_to_file(result::AbstractMatrix; file_name::AbstractString, header::Vector{Float64})
    table = Tables.table(result)
    header = string.(header)
    CSV.write(file_name, table, header=header)

    return nothing
end
