# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using DataFrames
using CSV
using CairoMakie
using Base

function plot_density(experiment_path)
    csv_data = CSV.File(joinpath(experiment_path, "mps_density.csv"))
    df = DataFrame(csv_data)
    # Convert DataFrame to 2D Array of Floats
    density = Matrix{Float64}(df)
        
    num_times, num_atoms = size(density)
    
    ts = [i for i in 0:num_times-1]
    xs = [i for i in 0:num_atoms-1]
    
    fig, ax, hm = heatmap(xs, ts, density',  colorrange = (0, 1), axis=(;title = "Density evolution: n_i(t)", xlabel = "atom index, i", ylabel="time step, t_j"))
    Colorbar(fig[:, end+1])  # equivalent
    fig_path = joinpath(experiment_path, "mps_density.png")
    save(fig_path, fig)

    fig, ax, plt = lines(xs, density[end, :],  axis=(;title = "Final density: n_i(T)", xlabel = "atom index, i", ylabel="density, n_i"))
    fig_path = joinpath(experiment_path, "final_density.png")
    save(fig_path, fig)
end

function plot_correlator(experiment_path)
    csv_data = CSV.File(joinpath(experiment_path, "correlator_zz.csv"))
    df = DataFrame(csv_data)
    # Convert DataFrame to 2D Array of Floats
    corr_zz = Matrix{Float64}(df)
        
    num_atoms = size(corr_zz, 1)
    
    xs = [i for i in 0:num_atoms-1]    

    fig, ax, hm = heatmap(xs, xs, corr_zz', axis=(;title = "Correlator: <Sz_i(T) Sz_j(T)>", xlabel = "atom index, i", ylabel="atom index, j"))

    Colorbar(fig[:, end+1], colorrange = (-.25, .25))  # equivalent
    fig_path = joinpath(experiment_path, "correlator_zz.png")
    save(fig_path, fig)
end

function plot_atoms(experiment_path)
    csv_data = CSV.File(joinpath(experiment_path, "atom_coordinates.csv"))
    df = DataFrame(csv_data)
    # Convert DataFrame to 2D Array of Floats
    coords = Matrix{Float64}(df)

    fig, ax, hm = scatter(coords[1, :], coords[2, :], axis=(;title = "Atom coordinates: (x, y)", xlabel = "atom coordinate, x", ylabel="atom coordinate, y"))

    fig_path = joinpath(experiment_path, "atom_coordinates.png")
    save(fig_path, fig)
end

function plot_bitstrings(experiment_path, max_vals=20)
    csv_data = CSV.File(joinpath(experiment_path, "mps_samples.csv"))
    df = DataFrame(csv_data)
    # Convert DataFrame to 2D Array of Floats
    bit_values = Matrix{Int}(df)
    bitstrings = [join(string.(row)) for row in eachrow(bit_values)]

    # Create a counter (dictionary) from the list
    counter = Dict{String, Int}()
    for str in bitstrings
        counter[str] = get(counter, str, 0) + 1
    end    

    # Create a list of counts corresponding to each string    
    ks = collect(keys(counter))
    vs = collect(values(counter))
    
    # convert bitstrings to integer representation
    # Sort by counts
    sorted_indxs = sortperm(vs, rev = true)
    num_vals = length(vs)

    fig, ax, plt = barplot(1:num_vals, vs[sorted_indxs], axis = (xticks = (1:num_vals, ks[sorted_indxs]),
                           title = "Sampled bitstrings: t=T", xlabelrotation=1), )
    fig_path = joinpath(experiment_path, "mps_samples.png")
    save(fig_path, fig)
end


function plot_all(experiment_path)
    plot_density(experiment_path)
    plot_correlator(experiment_path)
    plot_atoms(experiment_path)
    plot_bitstrings(experiment_path)
end
