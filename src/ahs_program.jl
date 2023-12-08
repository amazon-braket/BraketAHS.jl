# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using ITensors
using CSV, DataFrames
using Dates
using Base.Filesystem
using Missings
using Random
using JSON, JSON3

using Markdown
using InteractiveUtils

using Braket
using Braket: AtomArrangement, AtomArrangementItem, TimeSeries, DrivingField, AwsDevice, AnalogHamiltonianSimulation, discretize, AnalogHamiltonianSimulationQuantumTaskResult
using DataStructures, Statistics, Plots

begin
	a = 5.5e-6	
	register = AtomArrangement()
	push!(register, AtomArrangementItem((0.5, 0.5 + 1/√2) .* a))
	push!(register, AtomArrangementItem((0.5 + 1/√2, 0.5) .* a))
	push!(register, AtomArrangementItem((0.5 + 1/√2, -0.5) .* a))
	push!(register, AtomArrangementItem((0.5, -0.5 - 1/√2) .* a))
	push!(register, AtomArrangementItem((-0.5, -0.5 - 1/√2) .* a))
	push!(register, AtomArrangementItem((-0.5 -1/√2, -0.5) .* a))
	push!(register, AtomArrangementItem((-0.5 -1/√2, 0.5) .* a))
	push!(register, AtomArrangementItem((-0.5, 0.5 + 1/√2) .* a))
end


begin
	time_max                = 4e-6  # seconds
	time_ramp               = 1e-7  # seconds
	Ω_max                   = 6300000.0  # rad / sec
	Δ_start                 = -5 * Ω_max
	Δ_end                   = 5 * Ω_max
	
	Ω                       = TimeSeries()
	Ω[0.0]                  = 0.0
	Ω[time_ramp]            = Ω_max
	Ω[time_max - time_ramp] = Ω_max
	Ω[time_max]             = 0.0
	
	Δ                       = TimeSeries()
	Δ[0.0]                  = Δ_start
	Δ[time_ramp]            = Δ_start
	Δ[time_max - time_ramp] = Δ_end
	Δ[time_max]             = Δ_end
	
	ϕ           = TimeSeries()
	ϕ[0.0]      = 0.0
	ϕ[time_max] = 0.0
end

begin
	drive                   = DrivingField(Ω, ϕ, Δ)
	ahs_program             = AnalogHamiltonianSimulation(register, drive)
end


json_str = JSON3.write(ir(ahs_program))
json_obj = JSON.parse(json_str)

# Define the file path
file_path = "../examples/ahs_program.json"

# Write the JSON object to a file
open(file_path, "w") do file
    JSON.print(file, json_obj, 4)  # The '4' here is for pretty printing with an indent of 4 spaces
end


function get_atom_coordinates(lattice_sites, node_removal_prob)
    """ Generate atom positions for square lattice
    """
    # set random seed
    Random.seed!(1)
    # generate atom coords based on the lattice    
    N = length(lattice_sites)

    pruned_atom_coordinates = []
    filling = []
    for idx in 1:N
        if rand() > node_removal_prob
            push!(pruned_atom_coordinates, lattice_sites[idx])
            fill = 1
        else
            fill = 0
        end
        push!(filling, fill)
    end        
    return pruned_atom_coordinates, filling
end

