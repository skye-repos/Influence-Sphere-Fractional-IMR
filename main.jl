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
const N::Int = 2e4
const p::Float64 = 4.5 / (N - 1)
const num_realizations::Int = 25
const θ_list = collect(0.4:0.01:0.8)
const l_θ = length(θ_list)
const F₀_list = collect(0.10:0.05:0.40)
const l_F = length(F₀_list)

const graphs = make_graphs()
const plot_labels = ["F₀ = $(i)" for i ∈ F₀_list]

## This runs the actual simulations and averages over 25 realizations
@time flips, ffrac, exin0, exin1 = average_simulations(graphs);

exins = exin1 - exin0

exin0_rescaled = zeros(Float64, (l_θ, l_F))
exin1_rescaled = zeros(Float64, (l_θ, l_F))
exins_rescaled = zeros(Float64, (l_θ, l_F))
for i ∈ eachindex(F₀_list)
    local max0 = maximum(exin0[:, i])
    local max1 = maximum(exin1[:, i])
    local maxs = maximum(abs.(exins[:, i]))
    for j ∈ eachindex(θ_list)
        exin0_rescaled[j, i] = exin0[j, i] / max0
        exin1_rescaled[j, i] = exin1[j, i] / max1
        exins_rescaled[j, i] = exins[j, i] / maxs
    end
end

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
    title="Final fraction of nodes w/ Oᵥ = 0, N = $(N)",
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

results_plotter(abs.(exins), θ_list;
    ylabel="|⟨I₁⟩ - ⟨I₀⟩|",
    title="Absolute Comparative Instability of Opinions, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-exins-abs.png")

results_plotter(abs.(exins_rescaled), θ_list;
    ylabel="|⟨I₁⟩ - ⟨I₀⟩|",
    title="Absolute Rescaled Comparative Instability, N = $(N)",
    labels=plot_labels)
savefig(plt_dir * "$(N)-exins-rescaled-abs.png")

hmp_flips = heatmap(θ_list, F₀_list, flips',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="% flip-floppers")
savefig(hmp_flips, plt_dir * "$(N)-heatmap-flips.png")

hmp_ffrac = heatmap(θ_list, F₀_list, ffrac',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="% with final Oᵥ = 0")
savefig(hmp_ffrac, plt_dir * "$(N)-heatmap-ffrac.png")

hmp_exin0_rsc = heatmap(θ_list, F₀_list, exin0_rescaled',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="Rescaled ⟨I₀⟩")
savefig(hmp_exin0_rsc, plt_dir * "$(N)-heatmap-exin0-rescaled.png")

hmp_exin0 = heatmap(θ_list, F₀_list, exin0',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="⟨I₀⟩")
savefig(hmp_exin0, plt_dir * "$(N)-heatmap-exin0.png")

hmp_exin1_rsc = heatmap(θ_list, F₀_list, exin1_rescaled',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="Rescaled ⟨I₁⟩")
savefig(hmp_exin1_rsc, plt_dir * "$(N)-heatmap-exin1-rescaled.png")

hmp_exin1 = heatmap(θ_list, F₀_list, exin1',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="⟨I₁⟩")
savefig(hmp_exin1, plt_dir * "$(N)-heatmap-exin1.png")

hmp_exins_rsc = heatmap(θ_list, F₀_list, exins_rescaled',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="Rescaled ⟨I₁⟩ - ⟨I₀⟩")
savefig(hmp_exins_rsc, plt_dir * "$(N)-heatmap-exins-rescaled.png")

hmp_exins = heatmap(θ_list, F₀_list, exins',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="⟨I₁⟩ - ⟨I₀⟩")
savefig(hmp_exins, plt_dir * "$(N)-heatmap-exins.png")

hmp_exins_rsc_abs = heatmap(θ_list, F₀_list, abs.(exins_rescaled)',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="Rescaled |⟨I₁⟩ - ⟨I₀⟩|")
savefig(hmp_exins_rsc_abs, plt_dir * "$(N)-heatmap-exins-rescaled-abs.png")

hmp_exins_abs = heatmap(θ_list, F₀_list, abs.(exins)',
    c=:viridis,
    xlabel="θ", ylabel="F₀",
    title="|⟨I₁⟩ - ⟨I₀⟩|")
savefig(hmp_exins_abs, plt_dir * "$(N)-heatmap-exins-abs.png")
