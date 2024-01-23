folder_path = joinpath(splitpath(@__FILE__)[1:end-1])

# Define the file path
file_path = joinpath(folder_path, "mis.json")

# Define the folder to save plots and data
experiment_path = joinpath(folder_path, "mis_plots")

using TensorNetworkAHSimulator
using Braket
using JSON, JSON3


begin
    num_rows = 7 # odd number
    num_cols = num_rows # odd number
	a = 5.5e-6	
	register = AtomArrangement()
    for i in 1 : num_rows
        for j in 1 : num_cols
            push!(register, AtomArrangementItem((i * a, j * a)))
        end
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
