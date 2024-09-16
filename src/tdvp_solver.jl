using ITensorTDVP

function get_interaction_tdvp(Vij, N)
    H_int = OpSum()
    for j1 in 1:N, j2 in (j1+1):N
        if Vij[j1, j2] != 0
            # Recall the number operator n = 1/2+Sz
            # and H_int = Vij * ni * nj

            H_int += Vij[j1, j2], "Sz", j1, "Sz", j2
            H_int += Vij[j1, j2]/2, "I",  j1, "Sz", j2
            H_int += Vij[j1, j2]/2, "Sz", j1, "I",  j2
            H_int += Vij[j1, j2]/4, "I",  j1, "I",  j2
        end
    end
    return H_int
end

function get_drive_tdvp(protocol, i_τ, N)
    Ω_ts = protocol[:rabi_driving]
    Δ_glob_ts = protocol[:global_detuning]
    Δ_loc_ts = protocol[:local_detuning]
    pattern = protocol[:pattern]
    # Local Detuning: - Δ_loc_ts(t) * pattern[i] * n_i
    H_drive = OpSum()
    for j in 1:N
        H_drive += Ω_ts[i_τ], "Sx", j
        H_drive += -Δ_glob_ts[i_τ], "I", j
        H_drive += -Δ_glob_ts[i_τ], "Sz", j


        # H_drive += -Δ_loc_ts[i_τ]/2 * pattern[j], "I", j
        # H_drive += -Δ_loc_ts[i_τ]/2 * pattern[j], "Sz", j
    end
    return H_drive
end

function compute_MPS_evolution_tdvp(protocol, Vij, N, n_τ_steps, max_bond_dim, cutoff)
    τ = protocol[:τ]

    # s = siteinds("S=1/2", N; conserve_qns=false)
    s = siteinds("S=1/2", N; conserve_qns=false)

    # ψ = MPS(s, n -> "Dn")
    ψ = MPS(s, ["Dn" for _ in s])
    # ψ = random_mps(s, "↑"; linkdims=10)

    H_int = get_interaction_tdvp(Vij, N)
    for i_τ in 1:n_τ_steps
        # println(i_τ)
        H = H_int + get_drive_tdvp(protocol, i_τ, N)
        op = MPO(-im * τ * H, s)
        ψ = tdvp(op, 1, ψ; nsteps=1, maxdim=max_bond_dim, cutoff=cutoff, normalize=true)
    end

    return ψ
end