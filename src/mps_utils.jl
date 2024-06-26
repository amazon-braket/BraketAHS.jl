# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using ITensors
# using CSV, DataFrames
using Dates
using Base.Filesystem
using Missings
using Random
using JSON3


"""
Generate atom positions from the ahs program.
"""
function get_atom_coordinates(ahs_json)
    atom_coords = ahs_json["setup"]["ahs_register"]["sites"]
    filling = ahs_json["setup"]["ahs_register"]["filling"]
    atom_coords = [parse.(Float64, c) for c in atom_coords]
    return atom_coords, filling
end

function get_Vij(atom_coordinates, N::Int, interaction_R::Float64, C6::Float64)
    Vij = zeros(Float64, N, N)
    for (i, Ri) in enumerate(atom_coordinates)
        for (j, Rj) in enumerate(atom_coordinates)
            # NN approximation (truncate van der Waals tail)
            Rij = sqrt((Ri[1]-Rj[1])^2 + (Ri[2]-Rj[2])^2)
            # check if two atoms are close enough (but not zero)
            check_NN = 0 < Rij <= interaction_R
            check_NN && (Vij[i, j] = C6/abs(Rij)^6)
        end
    end
    return Vij
end

function piecewise_protocol(x, points, values)
    y = missing
    (x<points[1] || x>points[end]) && throw(ArgumentError("Input parameter x ($x) must be between min ($(points[1])) and max ($(points[end])) points values"))
    for (i, point) in enumerate(points[1:end-1])
        if x >= point && x <= points[i+1]
            k = (values[i+1]-values[i])/(points[i+1]-point)
            y = k*(x-point) + values[i]
        end
    end
    return y
end

function parse_protocol(ahs_program, n_τ_steps::Int)
    # Define piecewise functions (protocols)
    time_steps = collect(0:(n_τ_steps-1)) ./ n_τ_steps 
    time_points_Δ = ahs_program["hamiltonian"]["drivingFields"][1]["detuning"]["time_series"]["times"]
    values_Δ = ahs_program["hamiltonian"]["drivingFields"][1]["detuning"]["time_series"]["values"]

    time_points_Ω = ahs_program["hamiltonian"]["drivingFields"][1]["amplitude"]["time_series"]["times"]
    values_Ω = ahs_program["hamiltonian"]["drivingFields"][1]["amplitude"]["time_series"]["values"]

    if "localDetuning" ∉ keys(ahs_program["hamiltonian"]) || length(ahs_program["hamiltonian"]["localDetuning"]) == 0
        # Define an empty local detuning with zero pattern and zero values
        filling = ahs_program["setup"]["ahs_register"]["filling"]
        pattern = ["0" for _ in filling]
        time_points_local_detuning = ["0.0", time_points_Ω[end]]
        values_local_detuning = ["0.0", "0.0"]
    else
        pattern = ahs_program["hamiltonian"]["localDetuning"][1]["magnitude"]["pattern"]
        time_points_local_detuning = ahs_program["hamiltonian"]["localDetuning"][1]["magnitude"]["time_series"]["times"]
        values_local_detuning = ahs_program["hamiltonian"]["localDetuning"][1]["magnitude"]["time_series"]["values"]
    end

    # Convert strings to floats [parsing Braket AHS program]
    time_points_Δ = parse.(Float64, time_points_Δ)
    values_Δ = parse.(Float64, values_Δ)
    time_points_Ω = parse.(Float64, time_points_Ω)
    values_Ω = parse.(Float64, values_Ω)
    

    pattern = parse.(Float64, pattern)
    time_points_local_detuning = parse.(Float64, time_points_local_detuning)
    values_local_detuning = parse.(Float64, values_local_detuning)

    # Define piecewise protocols
    total_time = time_points_Δ[end]
    t_vals = [i/n_τ_steps*total_time for i in 1:n_τ_steps]
    # Global detuning time series
    Δ_glob_ts = [piecewise_protocol(t, time_points_Δ, values_Δ) for t in t_vals]
    # Local detuning time series
    Δ_loc_ts = [piecewise_protocol(t, time_points_local_detuning, values_local_detuning) for t in t_vals]
    # Rabi driving field
    Ω_ts = [piecewise_protocol(t, time_points_Ω, values_Ω) for t in t_vals]

    @assert length(time_steps) == n_τ_steps 
    @assert length(Ω_ts) == n_τ_steps 
    @assert length(Δ_loc_ts) == n_τ_steps 
    @assert length(Δ_glob_ts) == n_τ_steps 

    return (time_steps=time_steps,
            rabi_driving=Ω_ts,
            global_detuning=Δ_glob_ts,
            local_detuning=Δ_loc_ts,
            pattern=pattern,
            τ=total_time/n_τ_steps
            )
end

function compute_energies(samples, Vij::Matrix{Float64}, Δ_glob_ts, Δ_loc_ts, pattern)
    Δ_glob_end = Δ_glob_ts[end]
    Δ_loc_end = Δ_loc_ts[end]
    @info "Final global detuning (Δ_glob_end): $Δ_glob_end"
    @info "Final local detuning (Δ_loc_end): $Δ_loc_end"
    energies = map(eachrow(samples)) do sample
        energy_Vij = 0.5 * dot(sample, Vij, sample)
        energy_Δ_glob = -Δ_glob_end * sum(sample)
        # compute sum_i (pattern_i * n_i)
        energy_Δ_loc = -Δ_loc_end * dot(sample, pattern)
        return energy_Vij + energy_Δ_glob + energy_Δ_loc
    end
    return energies
end

"""
    get_trotterized_circuit_2d(sites, n_steps, N, Vij::Matrix{Float64})

Second order Trotterization circuit for a time-dependent Hamiltonian,
for a time step `n_steps` total time steps, on `N` total atoms.
`sites` defines the site indices in the MPS used to build the circuit.
Returns a `Vector{Vector{iTensor}}` list of gates at each time step.
"""
function get_trotterized_circuit_2d(sites, n_steps::Int, N::Int, Vij::Matrix{Float64}, protocol)
    τ = protocol[:τ]
    circuit = Vector{Vector{ITensor}}(undef, n_steps)
    for i_τ in 1:n_steps
        two_site_gates = ITensor[]
        ## Two-site terms: Vij*n_i*n_j
        for j1 in 1:N, j2 in (j1+1):N
            if Vij[j1, j2] != 0
                s1_op = (0.5*op("I", sites[j1]) + op("Sz", sites[j1]))
                s2_op = (0.5*op("I", sites[j2]) + op("Sz", sites[j2]))
                hj = Vij[j1, j2] * s1_op * s2_op
                Gj = exp(-im * τ / 2 * hj)
                push!(two_site_gates, Gj)
            end
        end
        Ω_ts = protocol[:rabi_driving]
        # Single-site terms: Ω_ts(t)*σX_i/2
        rabi_pulse_gates = map(1:N) do j
            hj_Ω = Ω_ts[i_τ] * op("Sx", sites[j])
            Gj = exp(-im * τ / 2 * hj_Ω)
            return Gj
        end

        Δ_glob_ts = protocol[:global_detuning]
        detuning_ops = [0.5*op("I", s) + op("Sz", s) for s in sites]
        # Global Detuning: - Δ_glob_ts(t)*n_i
        global_detuning_gates = map(detuning_ops) do op
            hj_Δ_glob = - Δ_glob_ts[i_τ] * op
            Gj = exp(-im * τ / 2 * hj_Δ_glob)
            return Gj
        end
        Δ_loc_ts = protocol[:local_detuning]
        pattern = protocol[:pattern]
        # Local Detuning: - Δ_loc_ts(t) * pattern[i] * n_i
        local_detuning_gates = map(zip(pattern, detuning_ops)) do (fill, op) 
            hj_Δ_loc = -Δ_loc_ts[i_τ] * fill * op
            Gj = exp(-im * τ / 2 * hj_Δ_loc)
            return Gj
        end
        all_gates = vcat(two_site_gates, rabi_pulse_gates, global_detuning_gates, local_detuning_gates)
        # Include gates in reverse order too: (N,N-1),(N-1,N-2),...
        append!(all_gates, reverse(all_gates))
        circuit[i_τ] = all_gates
    end
    return circuit
end

"""
    parse_ahs_program(parsed_args) -> (Vij, protocol, N)

Parse the AHS program stored in JSON format at `program_path` and prepare
experimental results directory at `experiment_path`.
Returns prepared experiment protocol `protocol` as a `NamedTuple` with keys
`(:time_steps, :rabi_driving, :global_detuning, :local_detuning, :pattern)`.
`Vij` is the generated inter-atomic potential, and `N` is the total number of atoms.
"""
function parse_ahs_program(ahs_json, args::Dict{String, Any})
    program_path = args["program-path"]
    experiment_path = args["experiment-path"]
    n_τ_steps = args["n-tau-steps"]
    interaction_R = args["interaction-radius"]
    C6 = args["C6"]

    # Read atom coords and fillings
    atom_coordinates, filling = get_atom_coordinates(ahs_json)
    N = length(atom_coordinates)

    if !isdir(experiment_path)
        mkdir(experiment_path)
        @info "Directory '$experiment_path' created."
    else
        @info "Directory '$experiment_path' already exists."
    end
    
    # Serialize the data to a JSON-formatted string, write the JSON string to the file
    json_str = JSON3.write(ahs_json)
    open(joinpath(experiment_path, "ahs_program.json"), "w") do file
        write(file, json_str)
    end

    # # Write filling data to CSV file
    # CSV.write(joinpath(experiment_path, "filling.csv"), DataFrame(filling', :auto))

    # array_2d_coords = hcat([t[1] for t in atom_coordinates], [t[2] for t in atom_coordinates])
    # CSV.write(joinpath(experiment_path, "atom_coordinates.csv"), DataFrame(array_2d_coords', :auto))
    # @debug "Atoms in atom array: $atom_coordinates"

    Vij = get_Vij(atom_coordinates, N, interaction_R, C6)
    protocol = parse_protocol(ahs_json, n_τ_steps)
    return Vij, protocol, N
end

"""
    compute_MPS_evolution(ψ::MPS,
                          circuit::Vector{Vector{ITensor}},
                          max_bond_dim::Int,
                          cutoff::Float64;
                          compute_truncation_error::Bool=false)

Evolve the MPS `ψ` according to the Trotterized circuit `circuit`.
`max_bond_dim` controls the maximum bond dimension
of the MPS, and `cutoff` controls the SVD truncation at each evolution step.
At each step, the per-atom `<Sz>` observable is computed and stored.
If `compute_truncation_error` is `true`, the error due to truncation is
computed at each step. The default is `false`.

!!! warning
    Computing the truncation error at every step is very computationally
    expensive.

Returns the `<Sz>` value on each atom at each time step, and the
truncation error at each time step
(`0.` at all times if `compute_truncation_error` is `false`),
and `ψ` after evolution.
"""
function compute_MPS_evolution(ψ::MPS, circuit::Vector{Vector{ITensor}}, max_bond_dim::Int, cutoff::Float64; compute_truncation_error::Bool=false)
    n_τ_steps = length(circuit)
    N = length(ψ)
    meas_array = Matrix{Float64}(undef, length(ψ), n_τ_steps)
    err_array = zeros(Float64, n_τ_steps)
    @info "Applying Trotter gates"
    for i_τ in 1:n_τ_steps
        # Save density expectation numbers
        Sz = expect(ψ, "Sz"; sites=1:N)
        meas_array[:, i_τ] = 0.5 .+ Sz
        # Evolution step: psi update
        ψ_prev = ψ
        ψ = ITensors.apply(circuit[i_τ], ψ; cutoff=cutoff, maxdim=max_bond_dim)
        normalize!(ψ)
        if compute_truncation_error
            ψ_exact = ITensors.apply(circuit[i_τ], ψ_prev; cutoff=1e-16)
            err_array[i_τ] = 1.0 - abs(inner(ψ_exact, ψ))
        end            
        @info "Step: $i_τ, current MPS bond dimension is $(maxlinkdim(ψ))"         
    end
    return meas_array, err_array, ψ    
end

function run(ahs_json, args)

    experiment_path = args["experiment-path"]
    n_τ_steps = args["n-tau-steps"]
    C6 = args["C6"]
    interaction_R = args["interaction-radius"]
    n_shots = args["shots"]    
    Vij, protocol, N = parse_ahs_program(ahs_json, args)

    @info "Preparing initial ψ MPS"
    s = siteinds("S=1/2", N; conserve_qns=false)

    # Initialize ψ to be a product state: Down state (Ground state)
    ψ = MPS(s, n -> "Dn")
    @info "Generating Trotterized circuit"
    circuit = get_trotterized_circuit_2d(s, n_τ_steps, N, Vij, protocol)

    max_bond_dim = args["max-bond-dim"]
    cutoff = args["cutoff"]
    compute_truncation_error = args["compute-truncation-error"]

    @info "Starting MPS evolution"
    res = @timed begin
        density, err_array, ψ = compute_MPS_evolution(ψ, circuit, max_bond_dim, cutoff, compute_truncation_error=compute_truncation_error)
    end
    # summary_array = ["time: $(res.time)",
    #                 "n_atoms: $N",
    #                 "Trotter steps: $n_τ_steps",
    #                 "interaction_R: $interaction_R",
    #                 "MPS cutoff: $cutoff",
    #                 "max_bond_dim: $max_bond_dim",
    #                 "total truncation err: $(sum(err_array))"]
    # summary_string = join(summary_array, ", ")
    # write(joinpath(experiment_path, "summary.txt"), summary_string)

    @info "Simulation complete. Elapsed time and memory used: $(res.time)."
    @info "Number of atoms: $N, MPS cutoff: $cutoff, max_bond_dim: $max_bond_dim, Trotter steps: $n_τ_steps"

    # Bitstring samples
    @info "Sampling from final MPS state"
    samples = Matrix{Int}(undef, N, n_shots)
    for shot in 1:n_shots
        sample_i = sample!(ψ) # Sampling bitstrings from a final psi(T)
        # iTensor MPS sample outputs values [1, 2] for 2-level system
        # Converting [1,2] -> [0,1]
        @views samples[:, shot] = [(2 - x) for x in sample_i]
    end

    return samples

    # # Correlation matrix 
    # correlator_zz = []

    # if args["compute-correlators"]
    #     @info "Evaluating correlation function ..."
    #     correlator_zz = 4 .* correlation_matrix(ψ, "Sz", "Sz") # renormalize to [-1, 1] range
    # end    

    # # Energies at t=T
    # energies = []
    
    # if args["compute-energies"]
    #     @info "Evaluating energies at t=T ..."

    #     Δ_glob_ts = protocol[:global_detuning]
    #     Δ_loc_ts = protocol[:local_detuning]
    #     pattern = protocol[:pattern]
    
    #     energies = compute_energies(samples', Vij, Δ_glob_ts, Δ_loc_ts, pattern)
    # end
    

    # results = Dict(
    #     "samples" => samples,
    #     "density" => density,
    #     "summary" => summary_array
    # )

    # if args["compute-energies"]
    #     results["energies"] = energies
    # end

    # if args["compute-correlators"]
    #     results["correlator_zz"] = correlator_zz
    # end    

    # return results
end

function run_batch(ahs_jsons, args; max_parallel=-1)
    # println("in run_batch")
    n_tasks = length(ahs_jsons)
    println(n_tasks)
    todo_tasks_ch = Channel{Int}(ch->foreach(ix->put!(ch, ix), 1:n_tasks), n_tasks)
    max_parallel_threads = max_parallel > 0 ? max_parallel : 32
    n_task_threads = min(max_parallel_threads, n_tasks)

    results = Vector{}(undef, n_tasks)
    function process_work()
        while isready(todo_tasks_ch)
            println("in while")
            my_ix = -1
            # need to lock the channel as it may become empty
            # and "unready" in between the while-loop call
            # and the call to take!
            lock(todo_tasks_ch) do
                my_ix = isready(todo_tasks_ch) ? take!(todo_tasks_ch) : -1
            end
            # if my_ix is still -1, the channel is empty and
            # there's no more work to do
            my_ix == -1 && break
            ahs_json  = ahs_jsons[my_ix]
            results[my_ix] = run(ahs_json, args)
        end
        return
    end
    tasks = Vector{}(undef, n_task_threads)
    @sync for worker in 1:n_task_threads
        println("worker = $worker")
        tasks[worker] = Threads.@spawn process_work()
    end
    # tasks don't return anything so we can wait rather than fetch
    wait.(tasks)
    # check to ensure all the results were in fact populated
    for r_ix in 1:n_tasks
        @assert isassigned(results, r_ix)
    end
    return results

end 


"""
    save_results(experiment_path::String,
                 ψ::MPS,
                 meas_array::Matrix{Float64},
                 n_shots::Int,
                 Vij::Matrix{Float64},
                 Δ_glob_ts)

Compute final energies and correlation functions from the MPS simulation, and
write them in `experiment_path`.

`experiment_path` is the path to the folder containing all data for this experiment.
`ψ` is the final MPS after time evolution. `meas_array` is the per-atom `Sz` measurement
computed at each time step in the simulation, returned from `compute_MPS_evolution`.
`n_shots` is the number of samples (shots) to compute from the final MPS `ψ`. `Vij` is
the interatomic potential generated by `get_Vij`. `Δ_glob_ts` is the global detuning function
time series.
"""

# function save_results(results, experiment_path::String)
#     @info "Saving results"
#     # # Samples
#     samples_file_path = joinpath(experiment_path, "mps_samples.csv") 
#     CSV.write(samples_file_path, DataFrame(results["samples"]', :auto), writeheader=false)
#     @info "Samples are saved to $samples_file_path"

#     # Energies
#     if haskey(results, "energies")
#         path = joinpath(experiment_path, "energies.csv")
#         CSV.write(path, DataFrame(results["energies"]', :auto), writeheader=false)
#         @info "Energies are saved to $path"
#     end

#     # Correlators
#     if haskey(results, "correlator_zz")
#         path = joinpath(experiment_path, "correlator_zz.csv")
#         CSV.write(path, DataFrame(real.(results["correlator_zz"]'), :auto), writeheader=false)
#         @info "Correlators are saved to $path"
#     end

#     # Densities
#     if haskey(results, "density")
#         path = joinpath(experiment_path, "mps_density.csv")
#         CSV.write(path, DataFrame(results["density"]', :auto), writeheader=false)
#         @info "Density evolution is saved to $path"
#     end
        
# end
