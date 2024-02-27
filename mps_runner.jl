# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using ArgParse
using ITensors
using CSV, DataFrames
using Dates
using Missings
using Random
using Logging
using JSON

include("src/mps_utils.jl")
include("src/plotter.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--program-path"
            help = "the path to the AHS program JSON file"
            arg_type = String
            required = true
        "--interaction-radius"
            help = "the interaction radius in meters"
            arg_type = Float64
            default = 7e-6
        "--experiment-path"
            help = "the directory in which to store all experiment data"
            arg_type = String
            default = joinpath("examples", "experiment_braket")
        "--cutoff"
            help = "cutoff for SVD values in MPS evolution"
            arg_type = Float64
            default = 1e-7
        "--shots"
            help = "number of shots for sampling"
            arg_type = Int
            default = 1000
        "--max-bond-dim"
            help = "maximum bond dimension for MPS"
            arg_type = Int
            default = 16
        "--compute-truncation-error"
            help = "whether to compute the error induced by truncation at each step (computationally expensive)"
            action = :store_true # default without this flag is false
        "--tau"
            help = "time evolution step size in seconds"
            arg_type = Float64
            default = 0.01e-6
        "--n-tau-steps"
            help = "number of time evolution steps to simulate"
            arg_type = Int
            default = 400
        "--C6"
            help = "C6 constant for van der Waals interaction between atoms in Rydberg state (Hz*m^6)"
            arg_type = Float64
            default = 5.42e-24
        "--compute-correlators"
            help = "Compute ZZ correlators at the end of the evolution (t=T)"
            action = :store_true
        "--compute-energies"
            help = "Compute energies from samples at the end of the evolution (t=T)"
            action = :store_true            
        "--generate-plots"
            help = "Generate plots after experiment is finished"
            action = :store_true
    end
    return parse_args(s)
end

function run(ahs_json, parsed_args)

    experiment_path = parsed_args["experiment-path"]
    τ = parsed_args["tau"]
    n_τ_steps = parsed_args["n-tau-steps"]
    C6 = parsed_args["C6"]
    interaction_R = parsed_args["interaction-radius"]
    n_shots = parsed_args["shots"]    
    Vij, protocol, N = parse_ahs_program(ahs_json, parsed_args)

    @info "Preparing initial ψ MPS"
    s = siteinds("S=1/2", N; conserve_qns=false)

    # Initialize ψ to be a product state: Down state (Ground state)
    ψ = MPS(s, n -> "Dn")
    @info "Generating Trotterized circuit"
    circuit = get_trotterized_circuit_2d(s, τ, n_τ_steps, N, Vij, protocol)

    max_bond_dim = parsed_args["max-bond-dim"]
    cutoff = parsed_args["cutoff"]
    compute_truncation_error = parsed_args["compute-truncation-error"]

    @info "Starting MPS evolution"
    res = @timed begin
        density, err_array, ψ = compute_MPS_evolution(ψ, circuit, max_bond_dim, cutoff, compute_truncation_error=compute_truncation_error)
    end
    summary_array = ["time: $(res.time)",
                    "n_atoms: $N",
                    "Trotter steps: $n_τ_steps",
                    "interaction_R: $interaction_R",
                    "MPS cutoff: $cutoff",
                    "max_bond_dim: $max_bond_dim",
                    "total truncation err: $(sum(err_array))"]
    summary_string = join(summary_array, ", ")
    write(joinpath(experiment_path, "summary.txt"), summary_string)

    @info "Simulation complete. Elapsed time and memory used: $(res.time)."
    @info "Number of atoms: $N, MPS cutoff: $cutoff, max_bond_dim: $max_bond_dim, Trotter steps: $n_τ_steps"

    # Bitstring samples
    @info "Sampling from final MPS state"
    samples = Matrix{Int}(undef, N, n_shots)
    for shot in 1:n_shots
        sample_i = sample!(ψ) # Sampling bitstrings from a final psi(T)
        # iTensor MPS sample outputs values [1, 2] for 2-level system
        # Converting [1,2] -> [0,1]
        samples[:, shot] = [(2 - val) for val in sample_i]
    end

    # Correlation matrix 
    correlator_zz = []

    if parsed_args["compute-correlators"]
        @info "Evaluating correlation function ..."
        correlator_zz = 4 .* correlation_matrix(ψ, "Sz", "Sz") # renormalize to [-1, 1] range
    end    

    # Energies at t=T
    energies = []
    
    if parsed_args["compute-energies"]
        @info "Evaluating energies at t=T ..."

        Δ_glob_ts = protocol[:global_detuning]
        Δ_loc_ts = protocol[:local_detuning]
        pattern = protocol[:pattern]
    
        energies = compute_energies(samples', Vij, Δ_glob_ts, Δ_loc_ts, pattern)
    end    

    results = Dict(
        "samples" => samples,
        "density" => density,
        "summary" => summary_array
    )

    if parsed_args["compute-energies"]
        results["energies"] = energies
    end

    if parsed_args["compute-correlators"]
        results["correlator_zz"] = correlator_zz
    end    

    return results
end

args = parse_commandline()
@info "Parsed command line arguments:"
for (k,v) in args
    @info "\t$k: $v"
end
experiment_path = args["experiment-path"]
program_path    = args["program-path"]

@info "JSON file to read: $program_path"
ahs_json = JSON.parsefile(program_path)

results = run(ahs_json, args)

@info "Saving results"
save_results(results, experiment_path)

@info "Generating plots"
if parsed_args["generate-plots"]
    @info "Plotting results from $experiment_path"
    plot_all(experiment_path)
    @info "Plotting complete."
end