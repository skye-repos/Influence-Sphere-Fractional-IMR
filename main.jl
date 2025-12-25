using Distributed
addprocs(8)
@everywhere begin
    include("utility.jl")
    include("simulate.jl")
    include("influence-radius.jl")
    using SharedArrays
end

@everywhere function average_simulations(
    graphs;
    θ_list=θ_list,
    F₀_list=F₀_list)
    
    l_g = num_realizations
    l_θ = length(θ_list)
    l_F = length(F₀_list)

    flips = SharedArray{Float64}(l_θ, l_F)
    ffrac = SharedArray{Float64}(l_θ, l_F)
    smode = SharedArray{Float64}(l_θ, l_F)

    @sync @distributed for n ∈ 1:l_g
        g = graphs[n]
        for (i, θ) ∈ pairs(θ_list), (j, F₀) ∈ pairs(F₀_list)
            local result = fractional_IMR(g, θ=θ, F₀=F₀)
            local p_s = calculate_instability(g, θ, result[1])
            flips[i, j] += count_flips_fraction(result) / l_g
            ffrac[i, j] += final_fraction(result) / l_g
            smode[i, j] += instability_mode(p_s) / l_g
            @info "Finished θ = $(θ), F₀ = $(F₀), and realization #$(n)/$(l_g)"
        end
    end

    return flips, ffrac, smode
end

## Declare constants
const N::Int = 1e3
const p::Float64 = 4.5 / (N - 1)
const m::Int = 5
const N0::Int = 15

const num_realizations::Int = 5
const type::String = "ER"

const θ_list = collect(0.05:0.01:0.95)
const F₀_list = collect(0.05:0.05:0.50)

const graphs = make_graphs(graph_dir)

const plot_labels = ["F₀ = $(i)" for i ∈ F₀_list]

## This runs the actual simulations and averages over 25 realizations
@time flips, ffrac, smode = average_simulations(graphs)

## Writing to CSV
CSV_dir = "./results/CSV/$(type)/N=$(N)/"
if !ispath(CSV_dir)
    mkpath(CSV_dir)
end

flips_CSV_file = CSV_dir * "flips " * join(string.(F₀_list), ", ") * ".csv"
write_to_file(flips',
    file_name=flips_CSV_file,
    header=θ_list)

ffrac_CSV_file = CSV_dir * "ffrac " * join(string.(F₀_list), ", ") * ".csv"
write_to_file(ffrac',
    file_name=ffrac_CSV_file,
    header=θ_list)

smode_CSV_file = CSV_dir * "smode " * join(string.(F₀_list), ", ") * ".csv"
write_to_file(smode',
    file_name=smode_CSV_file,
    header=θ_list)

## Plotting Everything
plot_dir = "./results/plots/$(type)/N=$(N)/";
if !ispath(plot_dir)
    mkpath(plot_dir)
end

smode_plot = results_plotter(smode;
    ylabel="% of nodes w/ max instability",
    title="N = $(N), $(type), Fractional-IMR m1/(m0+m1)")

savefig(plot_dir * "comparison-ismode.png")

flips_plot = results_plotter(flips;
    ylabel="% of nodes flipping > once",
    title="N = $(N), $(type), Fractional-IMR m1/(m0+m1)")

savefig(plot_dir * "comparison-#flips.png")

ffrac_plot = results_plotter(ffrac;
    ylabel="fraction of nodes with final state 0",
    title="N = $(N), $(type), Fractional-IMR m1/(m0+m1)")

savefig(plot_dir * "comparison-%ffrac.png")
