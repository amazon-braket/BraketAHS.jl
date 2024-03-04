using BraketAHS: run, save_results
using Test
using JSON, JSON3

@testset "BraketAHS.jl" begin
    # Write your tests here.

    ahs_json = JSON.parsefile("test/test_program.json")

    args = Dict(
        "experiment-path" => "test/",
        "program-path" => "",
        "interaction-radius" => 7e-6,
        "cutoff" => 1e-7,
        "shots" => 1000,
        "max-bond-dim" => 16,
        "compute-truncation-error" => false,
        "tau" => 0.01e-6,        
        "n-tau-steps" => 400,
        "C6" => 5.42e-24,
        "compute-correlators" => false,
        "compute-energies" => false,
        "generate-plots" => false
    )

    @show ahs_json

    result = run(ahs_json, args)

    mps_density = result["density"][:, end]
    expected_density = [0.49558400861051327,0.4955839949350386,0.49558370707943783,0.49558323080865946]
    @test mps_density == expected_density
    # add apprx cmoparison isapprox(array1, array2; atol=tolerance)
end
