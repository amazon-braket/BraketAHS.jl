using BraketAHS: run, save_results
using Test
using JSON3

@testset "BraketAHS.jl" begin
    # Write your tests here.

    json_str = """{"setup":{"ahs_register":{"sites":[["0.0","0.0"],["5.499999999999999856993733449161965e-6","0.0"],["0.0","5.499999999999999856993733449161965e-6"],["5.499999999999999856993733449161965e-6","5.499999999999999856993733449161965e-6"]],"filling":[1,1,1,1]}},"braketSchemaHeader":{"name":"braket.ir.ahs.program","version":"1"},"hamiltonian":{"shiftingFields":[{"magnitude":{"pattern":["0.0","0.0","0.0","0.0"],"time_series":{"values":["0.0","0.0"],"times":["0.0","3.999999999999999818992447303545035e-6"]}}}],"drivingFields":[{"phase":{"pattern":"uniform","time_series":{"values":["0.0","0.0"],"times":["0.0","3.999999999999999818992447303545035e-6"]}},"detuning":{"pattern":"uniform","time_series":{"values":["-3.15e7","-3.15e7","3.15e7","3.15e7"],"times":["0.0","9.999999999999999547481118258862587e-8","3.899999999999999929396754527743951e-6","3.999999999999999818992447303545035e-6"]}},"amplitude":{"pattern":"uniform","time_series":{"values":["0.0","6.3e6","6.3e6","0.0"],"times":["0.0","9.999999999999999547481118258862587e-8","3.899999999999999929396754527743951e-6","3.999999999999999818992447303545035e-6"]}}}]}}
    """
    ahs_json = JSON3.read(json_str)

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

    result = run(ahs_json, args)

    mps_density = result["density"][:, end]
    expected_density = [0.49558400861051327,0.4955839949350386,0.49558370707943783,0.49558323080865946]
    @test isapprox(mps_density, expected_density; atol=1e-7)
    # add apprx cmoparison isapprox(array1, array2; atol=tolerance)
end
