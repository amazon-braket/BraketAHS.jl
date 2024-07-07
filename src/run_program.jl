using DelimitedFiles

function run_program(
    program_path::String;
    interaction_radius::Float64 = 12e-6,
    cutoff::Float64 = 1e-7,
    shots::Int64 = 1000,
    max_bond_dim::Int64 = 4,
    n_tau_steps::Int64 = 80,
    C6::Float64 = 5.42e-24,
    )
    
    parsed_args = Dict(
        "program-path" => program_path,
        "interaction-radius" => interaction_radius,
        "cutoff" => cutoff,
        "shots" => shots,
        "max-bond-dim" => max_bond_dim,
        "n-tau-steps" => n_tau_steps,
        "C6" => C6,
    )
    
    @info "Parsed input arguments"
    for (k,v) in parsed_args
        @info "\t$k: $v"
    end

    @info "JSON file to read: $program_path"
    ahs_json = JSON3.read(read(program_path, String), Dict{String, Any})

    # parse inputs
    n_τ_steps = parsed_args["n-tau-steps"]
    C6 = parsed_args["C6"]
    interaction_R = parsed_args["interaction-radius"]
    Vij, protocol, N = parse_ahs_program(ahs_json, parsed_args)

    @info "Preparing initial ψ MPS"
    s = siteinds("S=1/2", N; conserve_qns=false)
    
    # Initialize ψ to be a product state: Down state (Ground state)
    ψ = MPS(s, n -> "Dn")
    @info "Generating Trotterized circuit"
    circuit = get_trotterized_circuit_2d(s, n_τ_steps, N, Vij, protocol)
    
    max_bond_dim = parsed_args["max-bond-dim"]
    cutoff = parsed_args["cutoff"]

    @info "Starting MPS evolution"
    res = @timed begin
        density, err_array, ψ = compute_MPS_evolution(ψ, circuit, max_bond_dim, cutoff, compute_truncation_error=false)
    end
    
    # Bitstring samples
    @info "Sampling from final MPS state"
    samples = Matrix{Int}(undef, N, shots)
    for shot in 1:shots
        sample_i = sample!(ψ) # Sampling bitstrings from a final psi(T)
        # iTensor MPS sample outputs values [1, 2] for 2-level system
        # Converting [1,2] -> [0,1]
        @views samples[:, shot] = [(2 - x) for x in sample_i]
    end

    experiment_path = joinpath(split(program_path,"/")[1:end-1])
    open(joinpath(experiment_path,"mps_samples.txt"), "w") do io
        writedlm(io, samples)
    end
end


# function run_program(
#     program_path::String;
#     interaction_radius::Float64 = 12e-6,
#     experiment_path::String = joinpath(split(program_path,"/")[1:end-1]),
#     cutoff::Float64 = 1e-7,
#     shots::Int64 = 1000,
#     max_bond_dim::Int64 = 4,
#     compute_truncation_error::Bool = false,
#     n_tau_steps::Int64 = 80,
#     C6::Float64 = 5.42e-24,
#     dim::Int64 = 2,
#     if_compute_correlators::Bool = false,
#     if_compute_energies::Bool = false,
#     )
    
#     # experiment_path = "/$(experiment_path)"
#     parsed_args = Dict(
#         "program-path" => program_path,
#         "interaction-radius" => interaction_radius,
#         "experiment-path" => experiment_path,
#         "cutoff" => cutoff,
#         "shots" => shots,
#         "max-bond-dim" => max_bond_dim,
#         "compute-truncation-error" => compute_truncation_error,
#         "n-tau-steps" => n_tau_steps,
#         "C6" => C6,
#         "dim" => dim,
#         "compute-correlators" => if_compute_correlators,
#         "compute-energies" => if_compute_energies,
#     )
    
#     @info "Parsed input arguments"
#     for (k,v) in parsed_args
#         @info "\t$k: $v"
#     end

#     @info "JSON file to read: $program_path"
#     ahs_json = JSON3.read(read(program_path, String), Dict{String, Any})

#     # parse inputs
#     experiment_path = parsed_args["experiment-path"]
#     n_τ_steps = parsed_args["n-tau-steps"]
#     C6 = parsed_args["C6"]
#     interaction_R = parsed_args["interaction-radius"]
#     Vij, protocol, N = parse_ahs_program(ahs_json, parsed_args)

#     @info "Preparing initial ψ MPS"
#     s = siteinds("S=1/2", N; conserve_qns=false)
    
#     # Initialize ψ to be a product state: Down state (Ground state)
#     ψ = MPS(s, n -> "Dn")
#     @info "Generating Trotterized circuit"
#     circuit = get_trotterized_circuit_2d(s, n_τ_steps, N, Vij, protocol)
    
#     max_bond_dim = parsed_args["max-bond-dim"]
#     cutoff = parsed_args["cutoff"]
#     compute_truncation_error = parsed_args["compute-truncation-error"]

#     @info "Starting MPS evolution"
#     res = @timed begin
#         density, err_array, ψ = compute_MPS_evolution(ψ, circuit, max_bond_dim, cutoff, compute_truncation_error=compute_truncation_error)
#     end
#     summary_array = ["time: $(res.time)",
#                       "n_atoms: $N",
#                       "Trotter steps: $n_τ_steps",
#                       "interaction_R: $interaction_R",
#                       "MPS cutoff: $cutoff",
#                       "max_bond_dim: $max_bond_dim",
#                       "total truncation err: $(sum(err_array))"]
#     summary_string = join(summary_array, ", ")
#     write(joinpath(experiment_path, "summary.txt"), summary_string)
    
#     @info "Simulation complete. Elapsed time and memory used: $(res.time)."
#     @info "Number of atoms: $N, MPS cutoff: $cutoff, max_bond_dim: $max_bond_dim, Trotter steps: $n_τ_steps"
    
    
#     # Bitstring samples
#     @info "Sampling from final MPS state"
#     samples = Matrix{Int}(undef, N, shots)
#     for shot in 1:shots
#         sample_i = sample!(ψ) # Sampling bitstrings from a final psi(T)
#         # iTensor MPS sample outputs values [1, 2] for 2-level system
#         # Converting [1,2] -> [0,1]
#         @views samples[:, shot] = [(2 - x) for x in sample_i]
#     end

#     # JSON3.write(joinpath(experiment_path, "mps_samples.json"),samples)

#     open(joinpath(experiment_path,"mps_samples.txt"), "w") do io
#         writedlm(io, samples)
#     end
 

#     # results = Dict(
#     #     "samples" => samples,
#     #     "density" => density,
#     #     "summary" => summary_array
#     # )

#     # if parsed_args["compute-correlators"]
#     #     @info "Evaluating correlation function ..."
#     #     correlator_zz = []
#     #     correlator_zz = 4 .* correlation_matrix(ψ, "Sz", "Sz") # renormalize to [-1, 1] range
#     #     results["correlator_zz"] = correlator_zz
#     # end    

#     # if parsed_args["compute-energies"]
#     #     @info "Evaluating energies at t=T ..."
#     #     energies = []

#     #     Δ_glob_ts = protocol[:global_detuning]
#     #     Δ_loc_ts = protocol[:local_detuning]
#     #     pattern = protocol[:pattern]
    
#     #     energies = compute_energies(samples', Vij, Δ_glob_ts, Δ_loc_ts, pattern)
#     #     results["energies"] = energies
#     # end

#     # save_results(results, experiment_path)
# end
