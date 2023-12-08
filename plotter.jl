# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using ArgParse
include("src/plotter.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--experiment-path"
            help = "the directory in which to store all experiment data"
            arg_type = String
            default = joinpath("examples", "experiment_braket")
    end
    return parse_args(s)
end


parsed_args = parse_commandline()
@info "Parsed command line arguments:"
for (k,v) in parsed_args
    @info "\t$k: $v"
end
experiment_path = parsed_args["experiment-path"]

@info "Plotting results from $experiment_path"
plot_all(experiment_path)
@info "Plotting complete."