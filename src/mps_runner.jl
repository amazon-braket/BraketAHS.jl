# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using ArgParse
using ITensors
using CSV, DataFrames
using Dates
using Missings
using Random
using Logging
using JSON3

include("mps_utils.jl")
include("plotter.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--number-of-atoms"
            help = "number of atoms in the program"
            arg_type = Int
            default = 16
        "--interaction-radius"
            help = "the interaction radius in meters"
            arg_type = Float64
            default = 7e-6
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
            default = 100
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
            action = :store_false
        "--compute-energies"
            help = "Compute energies from samples at the end of the evolution (t=T)"
            action = :store_false
        "--generate-plots"
            help = "Generate plots after experiment is finished"
            action = :store_false
        "--compute-density"
            help = "Compute and store density at each timestep"
            action = :store_false
    end
    return parse_args(s)
end

args = parse_commandline()
        
@info "Parsed command line arguments:"
for (k,v) in args
    @info "\t$k: $v"
end

args["experiment-path"] = "adiabatic_prep/results_$(args["number-of-atoms"])"
args["program-path"] = "adiabatic_prep/N_$(args["number-of-atoms"]).json"

@info "JSON file to read: $(args["program-path"])"
ahs_json = JSON3.read(read(args["program-path"], String), Dict{String, Any})

results = run(ahs_json, args)

@info "Saving results"
save_results(results, args["experiment-path"])

@info "Generating plots"
if args["generate-plots"]
    @info "Plotting results from $(args["experiment-path"])"
    plot_all(args["experiment-path"])
    @info "Plotting complete."
end
