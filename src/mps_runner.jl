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
        "--program-path"
            help = "the path to the AHS program JSON file"
            arg_type = String
            default = joinpath(dirname(@__DIR__), "examples", "ahs_program.json")
        "--interaction-radius"
            help = "the interaction radius in meters"
            arg_type = Float64
            default = 7e-6
        "--experiment-path"
            help = "the directory in which to store all experiment data"
            arg_type = String
            default = joinpath(dirname(@__DIR__), "examples", "experiment_braket")
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


args = parse_commandline()
        
@info "Parsed command line arguments:"
for (k,v) in args
    @info "\t$k: $v"
end

experiment_path = args["experiment-path"]
program_path    = args["program-path"]

@info "JSON file to read: $program_path"
ahs_json = JSON3.read(read(program_path, String), Dict{String, Any})

results = run(ahs_json, args)

@info "Saving results"
save_results(results, experiment_path)

@info "Generating plots"
if args["generate-plots"]
    @info "Plotting results from $experiment_path"
    plot_all(experiment_path)
    @info "Plotting complete."
end
