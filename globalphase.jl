using ITensors
using Braket
using Braket: AtomArrangement, AtomArrangementItem, TimeSeries, DrivingField, Pattern

begin
    a = 5.5e-6    
    register = AtomArrangement()
    push!(register, AtomArrangementItem((0., 0.) .* a))
    push!(register, AtomArrangementItem((1., 0.) .* a))
    push!(register, AtomArrangementItem((0., 1.) .* a))
    push!(register, AtomArrangementItem((1., 1.) .* a))
end

begin
    time_max                = 4e-6  # seconds
    time_ramp               = 1e-7  # seconds
    Ω_max                   = 6300000.0  # rad / sec
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

    Δ_loc                   = TimeSeries()
    Δ_loc[0.0]              = 0.
    Δ_loc[time_max]         = 0.

    ϕ                       = TimeSeries()
    ϕ[0.0]                  = 0.0
    ϕ[time_max]             = 0.0
end

begin
    # Modify to include global laser phase terms
    Φ_global                = 0.0  # Global laser phase term
    
    drive = DrivingField(Ω, ϕ .+ Φ_global, Δ)  # Include the global phase term
    
    pt = Pattern([0. for i in 1:length(register)])
    shift = ShiftingField(Field(Δ_loc, pt))

    ahs_program = AnalogHamiltonianSimulation(register, [drive, shift])
end


