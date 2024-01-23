folder_path = joinpath(splitpath(@__FILE__)[1:end-1])

# Define the file path
file_path = joinpath(folder_path, "LGT_1D_program.json")

# Define the folder to save plots and data
experiment_path = joinpath(folder_path, "LGT_1D_plots")

using TensorNetworkAHSimulator
using Braket
using JSON, JSON3

begin
    num_atoms = 45 # 4n+1
    num_defects = 2 # 1 or 2
	a = 5.5e-6	
	register = AtomArrangement()
    for i in 1 : num_atoms
        push!(register, AtomArrangementItem((0, (i-1) * a)))
    end
end


begin
    # t_max = 4e-6
    t_max = 8e-6
    Ω_ts = [
        0.0,
        1.0000000000000002e-06,
        1.0632000000000003e-06,
        1.1988349780753035e-06,
        1.2620349780753035e-06,
        2.262034978075304e-06,
        2.325234978075304e-06,
        3.325234978075304e-06,
        t_max
    ]

    Ω_vs = [0.0, 0.0, 15800000.0, 15800000.0, 0.0, 0.0, 15800000.0, 15800000.0, 0]

    shift_ts = [0,
    1.0000000000000002e-06,
    1.2620349780753038e-06,
    2.262034978075304e-06,
    2.325234978075304e-06,
    3.325234978075304e-06,
    t_max]
    shift_vs = [0, 47123889.803846896, 47123889.803846896, 0.0, 0, 0, 0]
	
	Ω                       = TimeSeries()
    for (t, v) in zip(Ω_ts, Ω_vs)
        Ω[t] = v
    end
	
	Δ                       = TimeSeries()
	Δ[0.0]                  = 0
    Δ[t_max]            = 0
	
	ϕ                       = TimeSeries()
	ϕ[0.0]                  = 0.0
	ϕ[t_max]             = 0.0

    shift_time_series = TimeSeries()
    for (t, v) in zip(shift_ts, shift_vs)
        shift_time_series[t] = v
    end
    
    if num_defects == 1
        pattern = zeros(Int, num_atoms)
        pattern[2:2:end] .= 1
        pattern[Int((num_atoms+1)/2)] = 1
    elseif num_defects == 2
        pattern = zeros(Int, num_atoms)
        pattern[2:2:end] .= 1
        ind1 = floor(Int, num_atoms/3)
        ind2 = floor(Int, num_atoms-num_atoms/3)
        pattern[ind1] == 0 ? pattern[ind1] = 1 : pattern[ind1+1] = 1
        pattern[ind2] == 0 ? pattern[ind2] = 1 : pattern[ind2-1] = 1
    end

    shift = Field(shift_time_series, Pattern(pattern))
end

begin
	drive                   = DrivingField(Ω, ϕ, Δ)
    shift                   = ShiftingField(shift)
	ahs_program             = AnalogHamiltonianSimulation(register, [drive, shift])
end


json_str = JSON3.write(ir(ahs_program))
json_obj = JSON.parse(json_str)



# Write the JSON object to a file
open(file_path, "w") do file
    JSON.print(file, json_obj, 4)  # The '4' here is for pretty printing with an indent of 4 spaces
end

run_program(file_path; experiment_path=experiment_path)
