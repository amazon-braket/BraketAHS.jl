# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using DataFrames
using CSV
using CairoMakie
using CairoMakie: heatmap, scatter, lines, Colorbar
using Base


function plot_density(experiment_path)
    csv_data = CSV.File(joinpath(experiment_path, "density.csv"), header=false)
    df = DataFrame(csv_data)
    density = Matrix{Float64}(df)
        
    num_times, num_atoms = size(density)
    
    ts = [i for i in 0:num_times-1]
    xs = [i for i in 0:num_atoms-1]
    
    fig, ax, hm = CairoMakie.heatmap(xs, ts, density', colorrange = (0, 1), axis=(;title = "Density evolution: n_i(t)", xlabel = "atom index, i", ylabel="time step, t_j"))
    Colorbar(fig[:, end+1])  # equivalent
    fig_path = joinpath(experiment_path, "density.png")
    save(fig_path, fig)

    fig, ax, plt = lines(xs, density[end, :],  axis=(;title = "Final density: n_i(T)", xlabel = "atom index, i", ylabel="density, n_i"))
    fig_path = joinpath(experiment_path, "final_density.png")
    save(fig_path, fig)
end

function plot_correlator(experiment_path)
    csv_data = CSV.File(joinpath(experiment_path, "correlator_zz.csv"), header=false)
    df = DataFrame(csv_data)
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
    coords = Matrix{Float64}(df)

    fig, ax, hm = scatter(coords[1, :], coords[2, :], axis=(;title = "Atom coordinates: (x, y)", xlabel = "atom coordinate, x", ylabel="atom coordinate, y"))

    fig_path = joinpath(experiment_path, "atom_coordinates.png")
    save(fig_path, fig)
end

function plot_bitstrings(experiment_path, max_vals=20)
    csv_data = CSV.File(joinpath(experiment_path, "bitstrings.csv"), header=false)
    df = DataFrame(csv_data)
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

function calculate_staggered_magnetization(mags)
    Nx = Ny = Int(sqrt(size(mags, 2)))
    num_times = size(mags, 1)
    ts = [i for i in 0:num_times-1]

    staggered_magnetization = zeros(Float64, num_times)
    for j in 1:Nx
        for i in 1:Ny
            idx_phys = Ny*(j-1) + i
            m = ( (-1)^(j+i - 2) ) * mags[:, idx_phys]
            staggered_magnetization += m
        end
    end
    return staggered_magnetization, ts, Nx, Ny
end

function variance_staggered_magnetization(mags, flattened_corr)
    staggered_magnetization, ts, Nx, Ny = calculate_staggered_magnetization(mags)
    var_stag_mag = -(staggered_magnetization .^ 2)

    n_timesteps = size(flattened_corr, 1)
    n_atoms = Int(sqrt(size(flattened_corr, 2)))  # Since it's n_atoms * n_atoms columns
    
    # Reshape into 3D array: timesteps × atoms × atoms
    corr_zz_3d = Array{Float64, 3}(undef, n_timesteps, n_atoms, n_atoms)
    for t in 1:n_timesteps
        corr_zz_3d[t, :, :] = reshape(flattened_corr[t, :], (n_atoms, n_atoms))
    end

    for t in 1:n_timesteps
        sum_corr = 0.0
        for j in 1:Nx
            for i in 1:Ny
                for jp in 1:Nx
                    for ip in 1:Ny
                        idx_phys_1 = Ny*(j-1) + i
                        idx_phys_2 = Ny*(jp-1) + ip
                        m = ( (-1)^(j+i+jp+ip - 4) ) * corr_zz_3d[t, idx_phys_1, idx_phys_2]
                        sum_corr += m
                    end
                end
            end
        end
        var_stag_mag[t] += sum_corr
    end
    var_stag_mag ./= (Nx*Ny)
    return var_stag_mag, ts
end

function plot_variance_staggered_magnetization(experiment_path)
    # Get density data for staggered magnetization
    csv_data = CSV.File(joinpath(experiment_path, "density.csv"), header=false)
    df = DataFrame(csv_data)
    mags = Matrix{Float64}(df)

    # Get correlator data and reconstruct 3D array
    csv_data = CSV.File(joinpath(experiment_path, "correlator_zz_timeseries.csv"), header=false)
    df = DataFrame(csv_data)
    flattened_corr = Matrix{Float64}(df)

    var_stag_mag, ts = variance_staggered_magnetization(mags, flattened_corr)
    
    fig, ax, plt = lines(ts, var_stag_mag, axis=(;title = "Variance of staggered magnetization: <(Sz_i(t) - <Sz_i(t)>)^2>/Nx*Ny", xlabel = "time step, t_j", ylabel="variance of staggered magnetization"))
    fig_path = joinpath(experiment_path, "variance_staggered_magnetization.png")
    save(fig_path, fig)
end

function plot_staggered_magnetization(experiment_path)
    csv_data = CSV.File(joinpath(experiment_path, "density.csv"), header=false)
    df = DataFrame(csv_data)
    mags = Matrix{Float64}(df)
    staggered_magnetization, ts, Nx, Ny = calculate_staggered_magnetization(mags)

    fig, ax, plt = lines(ts, staggered_magnetization, axis=(;title = "Staggered magnetization: <Sz_i(T)>", xlabel = "time step, t_j", ylabel="staggered magnetization, <Sz_i(T)>"))
    fig_path = joinpath(experiment_path, "staggered_magnetization.png")
    save(fig_path, fig)
end

function plot_correlator_evolution_gif(experiment_path)
    # Get correlator data and reconstruct 3D array
    csv_data = CSV.File(joinpath(experiment_path, "correlator_zz_timeseries.csv"), header=false)
    df = DataFrame(csv_data)
    flattened_corr = Matrix{Float64}(df)
    n_timesteps = size(flattened_corr, 1)
    n_atoms = Int(sqrt(size(flattened_corr, 2)))  # Since it's n_atoms * n_atoms columns
    
    # Reshape into 3D array: timesteps × atoms × atoms
    corr_zz_3d = Array{Float64, 3}(undef, n_timesteps, n_atoms, n_atoms)
    for t in 1:n_timesteps
        corr_zz_3d[t, :, :] = reshape(flattened_corr[t, :], (n_atoms, n_atoms))
    end
    
    xs = [i for i in 0:n_atoms-1]
    
    fig = Figure(resolution=(400, 450))
    ax = Axis(fig[1, 1], title="Correlator Evolution", xlabel="atom index, i", ylabel="atom index, j")
    hm = heatmap!(ax, xs, xs, corr_zz_3d[1, :, :], colorrange=(-0.25, 0.25))
    Colorbar(fig[1, 2], hm)
    
    # Create animation
    framerate = 10
    timestamps = Observable(1)
    
    on(timestamps) do t
        hm[3] = corr_zz_3d[t, :, :]
        ax.title = "t = $(t-1)"
    end
    
    # Record
    record(fig, joinpath(experiment_path, "correlator_zz_evolution.mp4"), 1:n_timesteps;
           framerate=framerate) do t
        timestamps[] = t
    end
end

function plot_many_variance_staggered_magnetization(experiment_paths)
    fig = Figure()
    ax = Axis(fig[1, 1], 
              title="Variance of staggered magnetization", 
              xlabel="time step, t_j", 
              ylabel="variance of staggered magnetization")
    
    for experiment_path in experiment_paths
        println(experiment_path)
        csv_data = CSV.File(joinpath(experiment_path, "density.csv"), header=false)
        df = DataFrame(csv_data)
        mags = Matrix{Float64}(df)

        # Get correlator data and reconstruct 3D array
        csv_data = CSV.File(joinpath(experiment_path, "correlator_zz_timeseries.csv"), header=false)
        df = DataFrame(csv_data)
        flattened_corr = Matrix{Float64}(df)

        var_stag_mag, ts = variance_staggered_magnetization(mags, flattened_corr)

        # Extract the system size from the path for the legend
        system_size = split(basename(experiment_path), "_")[2]
        lines!(ax, ts, var_stag_mag, label="N = $system_size")
    end
    
    axislegend(ax)  # Add the legend
    
    # Save the figure
    filename = joinpath(dirname(first(experiment_paths)), "variance_staggered_magnetization_comparison.png")
    println(filename)
    save(filename, fig)
    
    return fig
end

function plot_all(experiment_path)
    plot_density(experiment_path)
    plot_correlator(experiment_path)
    plot_correlator_evolution_gif(experiment_path)
    plot_atoms(experiment_path)
    plot_bitstrings(experiment_path)
    plot_staggered_magnetization(experiment_path)
    plot_variance_staggered_magnetization(experiment_path)
end

# Add main entry point to handle command line arguments
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 1
        println("Usage: julia --project=. src/plotter.jl <experiment_path>")
        exit(1)
    end
    experiment_path = ARGS[1]
    plot_all(experiment_path)
end

