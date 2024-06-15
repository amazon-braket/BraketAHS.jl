function run_program(
    program_path::String;
    interaction_radius::Float64 = 7e-6,
    experiment_path::String = joinpath("examples", "experiment_braket"),
    cutoff::Float64 = 1e-7,
    shots::Int64 = 1000,
    max_bond_dim::Int64 = 16,
    compute_truncation_error::Bool = false,
    n_tau_steps::Int64 = 400,
    C6::Float64 = 5.42e-24,
    dim::Int64 = 2,
    generate_plots::Bool = true
    )
    
    parsed_args = Dict(
        "program-path" => program_path,
        "interaction-radius" => interaction_radius,
        "experiment-path" => experiment_path,
        "cutoff" => cutoff,
        "shots" => shots,
        "max-bond-dim" => max_bond_dim,
        "compute-truncation-error" => compute_truncation_error,
        "n-tau-steps" => n_tau_steps,
        "C6" => C6,
        "dim" => dim,
        "generate-plots" => generate_plots
    )
    
    @info "Parsed input arguments"
    for (k,v) in parsed_args
        @info "\t$k: $v"
    end

    # parse inputs
    experiment_path = parsed_args["experiment-path"]
    n_τ_steps = parsed_args["n-tau-steps"]
    C6 = parsed_args["C6"]
    interaction_R = parsed_args["interaction-radius"]
    Vij, protocol, N = parse_ahs_program(parsed_args)

    @info "Preparing initial ψ MPS"
    s = siteinds("S=1/2", N; conserve_qns=false)
    
    # Initialize ψ to be a product state: Down state (Ground state)
    ψ = MPS(s, n -> "Dn")
    @info "Generating Trotterized circuit"
    circuit = get_trotterized_circuit_2d(s, τ, n_τ_steps, N, Vij, protocol)
    
    max_bond_dim = parsed_args["max-bond-dim"]
    cutoff = parsed_args["cutoff"]
    compute_truncation_error = parsed_args["compute-truncation-error"]

    @info "Starting MPS evolution"
    res = @timed begin
        meas_array, err_array, ψ = compute_MPS_evolution(ψ, circuit, max_bond_dim, cutoff, compute_truncation_error=compute_truncation_error)
    end
    summary_array = ["time: $(res.time)",
                      "n_atoms: $N",
                      "Trotter steps: $n_τ_steps",
                      "interaction_R: $interaction_R",
                      "MPS cutoff: $cutoff",
                      "max_bond_dim: $max_bond_dim",
                      "total truncation err: $(sum(err_array))"]
    summary_string = join(summary_array, ", ")
    write(joinpath(experiment_path, "summary.txt"), summary_string)
    
    @info "Elapsed time and memory used: $(res.time)."
    @info "Number of atoms: $N, MPS cutoff: $cutoff, max_bond_dim: $max_bond_dim, Trotter steps: $n_τ_steps"
    @info "Computing final correlation functions and saving results."
    
    save_results(experiment_path, ψ, meas_array, parsed_args["shots"], Vij, protocol, N)
    @info "Simulation complete."
    
    if parsed_args["generate-plots"]
        @info "Plotting results from $experiment_path"
        plot_all(experiment_path)
        @info "Plotting complete."
    end        

end