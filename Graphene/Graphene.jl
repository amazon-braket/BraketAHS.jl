folder_path = joinpath(splitpath(@__FILE__)[1:end-1])

# Define the file path
file_path = joinpath(folder_path, "graphene_program.json")

# Define the folder to save plots and data
experiment_path = joinpath(folder_path, "graphene_plots")

using TensorNetworkAHSimulator
using Braket
using JSON, JSON3

begin
    coords = [[0.5, -2.59807621],
    [3.5, -0.866025405],
    [-2.5, -0.866025405],
    [0.5, 0.866025405],
    [3.5, 2.59807621],
    [-2.5, 2.59807621],
    [0.5, 4.330127015],
    [-1.0, -3.464101615],
    [2.0, -1.732050805],
    [-4.0, -1.732050805],
    [-1.0, -6.80803795e-16],
    [2.0, 1.732050805],
    [-4.0, 1.732050805],
    [-1.0, 3.464101615],
    [2.0, -3.464101615],
    [-1.0, -1.732050805],
    [2.0, -1.37152845e-15],
    [-4.0, -1.364436935e-15],
    [-1.0, 1.732050805],
    [2.0, 3.464101615],
    [0.5, -4.330127015],
    [3.5, -2.59807621],
    [-2.5, -2.59807621],
    [0.5, -0.866025405],
    [3.5, 0.866025405],
    [-2.5, 0.866025405],
    [0.5, 2.59807621],
    [-0.5, -2.59807621],
    [2.5, -0.866025405],
    [-3.5, -0.866025405],
    [-0.5, 0.866025405],
    [2.5, 2.59807621],
    [-3.5, 2.59807621],
    [-0.5, 4.330127015],
    [-2.0, -3.464101615],
    [1.0, -1.732050805],
    [4.0, -5.96955965e-16],
    [-2.0, -2.04214679e-15],
    [1.0, 1.732050805],
    [-2.0, 3.464101615],
    [1.0, -3.464101615],
    [4.0, -1.732050805],
    [-2.0, -1.732050805],
    [1.0, -2.031550585e-15],
    [4.0, 1.732050805],
    [-2.0, 1.732050805],
    [1.0, 3.464101615],
    [-0.5, -4.330127015],
    [2.5, -2.59807621],
    [-3.5, -2.59807621],
    [-0.5, -0.866025405],
    [2.5, 0.866025405],
    [-3.5, 0.866025405],
    [-0.5, 2.59807621]]

	a = 5.75760852731661e-06
	register = AtomArrangement()
    for coord in coords
        push!(register, AtomArrangementItem((coord[1] * a, coord[2] * a)))
    end
end


begin
	time_max                = 4e-6  # seconds
	time_ramp               = 5e-7  # seconds
	Ω_max                   = 15800000.0  # rad / sec
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
    
    shift_time_series = TimeSeries()
    shift_time_series[0.0]      = 0.0
	shift_time_series[time_max] = 0.0
end

begin
	drive                   = DrivingField(Ω, ϕ, Δ)
    shift                   = ShiftingField(Field(shift_time_series, Pattern(zeros(length(register)))))
	ahs_program             = AnalogHamiltonianSimulation(register, [drive, shift])
end

json_str = JSON3.write(ir(ahs_program))
json_obj = JSON.parse(json_str)

# Write the JSON object to a file
open(file_path, "w") do file
    JSON.print(file, json_obj, 4)  # The '4' here is for pretty printing with an indent of 4 spaces
end

# Set interaction_radius=a effectively solve an MIS problem
run_program(file_path; experiment_path=experiment_path)
