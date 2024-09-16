
include("../src/mps_utils.jl")
include("../src/tdvp_solver.jl")

program_path = "examples/ahs_program_7q.json"
interaction_radius = 0.0 # 8.2e-6
cutoff = 1e-7
shots = 1000
max_bond_dim = 16
n_τ_steps = 100
C6 = 5.42e-24

parsed_args = Dict(
    "program-path" => program_path,
    "interaction-radius" => interaction_radius,
    "cutoff" => cutoff,
    "shots" => shots,
    "max-bond-dim" => max_bond_dim,
    "n-tau-steps" => n_τ_steps,
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


###### tebd
@info "Preparing initial ψ MPS"
s = siteinds("S=1/2", N; conserve_qns=false)

# Initialize ψ to be a product state: Down state (Ground state)
ψ1 = MPS(s, n -> "Dn")
@info "Generating Trotterized circuit"
circuit = get_trotterized_circuit_2d(s, n_τ_steps, N, Vij, protocol)

max_bond_dim = parsed_args["max-bond-dim"]
cutoff = parsed_args["cutoff"]

@info "Starting MPS evolution"
_, _, ψ2 = compute_MPS_evolution(ψ1, circuit, max_bond_dim, cutoff, compute_truncation_error=false)
samples_tebd = hcat([sample!(ψ2) for _ in 1:shots]...)
samples_tebd[samples_tebd.==1] .= 0
samples_tebd[samples_tebd.==2] .= 1

###### tdvp
ψ3 = compute_MPS_evolution_tdvp(protocol, Vij, N, n_τ_steps, max_bond_dim, cutoff)
samples_tdvp = hcat([sample!(ψ3) for _ in 1:shots]...)
samples_tdvp[samples_tdvp.==1] .= 0
samples_tdvp[samples_tdvp.==2] .= 1

display(samples_tebd)
display(samples_tdvp)
display(1 .- sum(samples_tebd, dims=2)./shots)
display(1 .- sum(samples_tdvp, dims=2)./shots)
println()
display(ψ2.data)
println()
display(ψ3.data)
