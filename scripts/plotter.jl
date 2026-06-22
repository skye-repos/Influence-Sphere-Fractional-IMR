using Pkg
Pkg.activate(".")
using Glob
using Plots
using CSV, DataFrames

files = glob("*/*/*", "../results/", join = true)

function make_heatmap(file::AbstractString, name::AbstractString)
    df = CSV.read(file, DataFrame)
    theta_list = parse.(Float64, names(df)[2:end])
    F0_list = df[:, 1]
    data = Matrix(df[:, 2:end])

    return heatmap(theta_list, F0_list, data, c = :viridis, title = name, xlabel = "θ values", ylabel = "F0 values")
end

out = "../plots/"
for file in files
    name = splitext(basename(file))[1]
    title = "placeholder"
    if name == "ctime"
        title = "Convergence Time"
    elseif name == "flips"
        title = "Fraction of flip-floppers"
    elseif name == "ffrac"
        title = "Fraction of nodes with final O(v) = 0"
    elseif name == "exin0"
        title = "Expected Instability of nodes with initial O(v) = 0"
    elseif name == "exin1"
        title = "Expected Instability of nodes with initial O(v) = 1"
    elseif name == "exins"
        title = "⟨I1⟩ - ⟨I0⟩"
    elseif name == "exin0_rsc"
        title = "Rescaled ⟨I0⟩"
    elseif name == "exin1_rsc"
        title = "Rescaled ⟨I1⟩"
    elseif name == "exins_rsc"
        title = "Rescaled ⟨I1⟩ - ⟨I0⟩"
    elseif name == "brnch"
        title = "Expected value of branching factor"
    elseif name == "outdg"
        title = "⟨k_out⟩ of Influence DiGraph"
    end

    hmp = make_heatmap(file, title)
    root = splitpath(file)[3]
    out_dir = out * root * "/"
    if !ispath(out_dir)
        mkpath(out_dir)
    end
    out_file = out_dir * name * ".png"
    savefig(hmp, out_file)
end
