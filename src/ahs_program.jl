# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using ITensors
using CSV, DataFrames
using Dates
using Base.Filesystem
using Missings
using Random
using JSON3

using Markdown
using InteractiveUtils

using Braket
using Braket: AtomArrangement, AtomArrangementItem, TimeSeries, DrivingField, Pattern
using DataStructures, Statistics, Plots

begin
	a = 5.5e-6	
	register = AtomArrangement()
	push!(register, AtomArrangementItem((0., 0.) .* a))
	push!(register, AtomArrangementItem((1., 0.) .* a))
	push!(register, AtomArrangementItem((0., 1.) .* a))
	push!(register, AtomArrangementItem((1., 1.) .* a))
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

	Δ_loc                   = TimeSeries()
	Δ_loc[0.0]              = 0.
	Δ_loc[time_max]             	= 0.

	ϕ           = TimeSeries()
	ϕ[0.0]      = 0.0
	ϕ[time_max] = 0.0
end

begin
	drive = DrivingField(Ω, ϕ, Δ)

	pt = Pattern([0. for i in 1:length(register)])
	shift = ShiftingField(Field(Δ_loc, pt))

	ahs_program = AnalogHamiltonianSimulation(register, [drive, shift])
end


json_str = JSON3.write(ir(ahs_program))
json_obj = JSON3.read(json_str)

# Define the file path
file_path = "examples/ahs_program.json"

# Write the JSON object to a file
open(file_path, "w") do file
    JSON3.write(file, json_obj)  # The '4' here is for pretty printing with an indent of 4 spaces
end

