using Pkg
Pkg.activate(".")
using Distributed
addprocs(10)
@everywhere begin
    include("lib/utility.jl")
    include("simulate.jl")
    include("influence.jl")
    using SharedArrays
end

@everywhere """
    average_simulations(
        graphs::AbstractVector;
        θ_list=θ_list,
        F₀_list=F₀_list
        )

Average simulations over num_realizations version of an ER graph. 
"""
function average_simulations(
    graphs::AbstractVector;
    θ_list=θ_list,
    F₀_list=F₀_list
)

    l_g = num_realizations
    l_θ = length(θ_list)
    l_F = length(F₀_list)

    flips = SharedArray{Float64}(l_θ, l_F)
    ffrac = SharedArray{Float64}(l_θ, l_F)
    exin0 = SharedArray{Float64}(l_θ, l_F)
    exin1 = SharedArray{Float64}(l_θ, l_F)

    @sync @distributed for n ∈ 1:l_g
        g = graphs[n]
        for (i, θ) ∈ pairs(θ_list), (j, F₀) ∈ pairs(F₀_list)
            result = fractional_IMR(g, θ=θ, F₀=F₀)
            Iₜ = total_instability(g, θ, result[end])
            I₀, I₁ = opinion_instability(Iₜ, result[end])
            flips[i, j] += count_flips(result) / l_g
            ffrac[i, j] += count_ffrac(result) / l_g
            exin0[i, j] += expected_instability(I₀) / l_g
            exin1[i, j] += expected_instability(I₁) / l_g
            @info "Finished θ = $(θ), F₀ = $(F₀), and realization #$(n)/$(l_g)"
        end
    end

    return flips, ffrac, exin0, exin1
end

## Declare simulation parameters
const N::Int = 3.5e4
const p::Float64 = 4.5 / (N - 1)
const num_realizations::Int = 25
const θ_list = collect(0.1:0.01:0.9)
const l_θ = length(θ_list)
const F₀_list = collect(0.15:0.05:0.45)
const l_F = length(F₀_list)

const graphs = make_graphs()
const plot_labels = ["F₀ = $(i)" for i ∈ F₀_list]

## This runs the actual simulations and averages over 25 realizations
@time flips, ffrac, exin0, exin1 = average_simulations(graphs);

exin0_rescaled = zeros(Float64, (l_θ, l_F))
exin1_rescaled = zeros(Float64, (l_θ, l_F))
for i ∈ eachindex(F₀_list)
    local max0 = maximum(exin0[:, i])
    local max1 = maximum(exin1[:, i])
    for j ∈ eachindex(θ_list)
        exin0_rescaled[j, i] = exin0[j, i] / max0
        exin1_rescaled[j, i] = exin1[j, i] / max1
    end
end

exins = exin1 - exin0
abs_exins = abs.(exins)
exins_rescaled = exin1_rescaled - exin0_rescaled
abs_exins_rescaled = abs.(exins_rescaled)

## Plotting & Writing to CSV
CSV_dir = "./results/CSV/N=$(N)/"
plt_dir = "./results/plots/N=$(N)/"

if !ispath(CSV_dir)
    mkpath(CSV_dir)
end

if !ispath(plt_dir)
    mkpath(plt_dir)
end

write_to_file(flips';
    directory=CSV_dir,
    file_name="flips",
    row_names=F₀_list,
    col_names=θ_list
)

write_to_file(ffrac';
    directory=CSV_dir,
    file_name="ffrac",
    row_names=F₀_list,
    col_names=θ_list
)

write_to_file(exin0';
    directory=CSV_dir,
    file_name="exin0",
    row_names=F₀_list,
    col_names=θ_list
)

write_to_file(exin1';
    directory=CSV_dir,
    file_name="exin1",
    row_names=F₀_list,
    col_names=θ_list
)

write_to_file(exins';
    directory=CSV_dir,
    file_name="exins",
    row_names=F₀_list,
    col_names=θ_list
)

results_plotter(flips, θ_list;
    ylabel="% nodes",
    title="Fraction of nodes that flip-flop, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-flips.png")

results_plotter(ffrac, θ_list;
    ylabel="% nodes",
    title="Equilibrium fraction of Oᵥ = 0, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-ffrac.png")

results_plotter(exin0, θ_list;
    ylabel="⟨I₀⟩",
    title="Expected Instability of Oᵥ = 0, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-exin0.png")

results_plotter(exin0_rescaled, θ_list;
    ylabel="⟨I₀⟩",
    title="Rescaled Instability of Oᵥ = 0, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-exin0-rescaled.png")

results_plotter(exin1, θ_list;
    ylabel="⟨I₁⟩",
    title="Expected Instability of Oᵥ = 1, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-exin1.png")

results_plotter(exin1_rescaled, θ_list;
    ylabel="⟨I₁⟩",
    title="Rescaled Instability of Oᵥ = 1, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-exin1-rescaled.png")

results_plotter(exins, θ_list;
    ylabel="⟨I₁⟩ - ⟨I₀⟩",
    title="Comparative Instability of Opinions, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-exins.png")

results_plotter(abs_exins, θ_list;
    ylabel="|⟨I₁⟩ - ⟨I₀⟩|",
    title="Magnitude of Comparative Instability, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-abs-exins.png")

results_plotter(exins_rescaled, θ_list;
    ylabel="⟨I₁⟩ - ⟨I₀⟩",
    title="Rescaled Comparative Instability, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-exins-rescaled.png")

results_plotter(abs_exins_rescaled, θ_list;
    ylabel="|⟨I₁⟩ - ⟨I₀⟩|",
    title="Magnitude of Rescaled Comparative Instability, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-abs-exins-rescaled.png")
