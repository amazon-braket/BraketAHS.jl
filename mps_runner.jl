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
        "--dim"
            help = "Dimension of lattice -- 1D or 2D"
            arg_type = Int
            default = 2
        "--generate-plots"
            help = "Generate plots after experiment is finished"
            action = :store_true
    end
    return parse_args(s)
end



parsed_args = parse_commandline()
@info "Parsed command line arguments:"
for (k,v) in parsed_args
    @info "\t$k: $v"
end
experiment_path = parsed_args["experiment-path"]
τ = parsed_args["tau"]
n_τ_steps = parsed_args["n-tau-steps"]
C6 = parsed_args["C6"]
interaction_R = parsed_args["interaction-radius"]
Vij, protocol, N = parse_ahs_program(parsed_args)

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
    meas_array, err_array, ψ = compute_MPS_evolution(ψ, circuit, max_bond_dim, cutoff, compute_truncation_error=compute_truncation_error)
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

@info "Elapsed time and memory used: $(res.time)."
@info "Number of atoms: $N, MPS cutoff: $cutoff, max_bond_dim: $max_bond_dim, Trotter steps: $n_τ_steps"
@info "Computing final correlation functions and saving results."

save_results(experiment_path, ψ, meas_array, parsed_args["shots"], Vij, protocol)
@info "Simulation complete."

if parsed_args["generate-plots"]
    @info "Plotting results from $experiment_path"
    plot_all(experiment_path)
    @info "Plotting complete."
end    