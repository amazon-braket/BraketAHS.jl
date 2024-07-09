using DelimitedFiles

function run_program(
    program_path::String;
    interaction_radius::Float64 = 9.2e-6,
    cutoff::Float64 = 1e-7,
    shots::Int64 = 1000,
    max_bond_dim::Int64 = 4,
    n_tau_steps::Int64 = 40,
    C6::Float64 = 5.42e-24,
    )

    shots = 1000

    
    println("wtf julia")
    println([shots, n_tau_steps, C6, interaction_radius, max_bond_dim])
        
    
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
    # Vij, protocol, N = parse_ahs_program(ahs_json, parsed_args)
    
    # Read atom coords and fillings
    atom_coordinates, _ = get_atom_coordinates(ahs_json)
    N = length(atom_coordinates)

    Vij = get_Vij(atom_coordinates, N, interaction_R, C6)
    protocol = parse_protocol(ahs_json, n_τ_steps)    

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

    open("mps_samples.txt", "w") do io
        writedlm(io, samples)
    end    
end
