include("thermodynamic_properties.jl")

mutable struct Mixture
    atoms::Array{Atom,1}
    molecules::Array{Molecule,1}
    vibrational_distribution::VibrationalDistribution

    Zr::Array{Float64,1}
    Er::Array{Float64,1}
    Zv::Array{Float64,1}
    Ev::Array{Float64,1}

    c_rot::Array{Float64,1}
    c_vibrT::Array{Float64,1}
    c_vibrTv::Array{Float64,1}

    c_v::Float64
    c_p::Float64

    U::Float64
    h::Float64

    c_vibr_equilibrium::Array{Float64,1}
    c_v_equilibrium::Float64
    c_p_equilibrium::Float64
end

function create_mixture(atoms, molecules, vd)
    return Mixture(atoms, molecules, vd,
    zeros(Float64, length(molecules)), zeros(Float64, length(molecules)),  # Zr/Er
    zeros(Float64, length(molecules)), zeros(Float64, length(molecules)),  # Zv/Ev
    zeros(Float64, length(molecules)), zeros(Float64, length(molecules)), zeros(Float64, length(molecules)), # c_rot, c_vibrT, c_vibrTv
    0.0, 0.0, 0.0, 0.0, # c_v, c_p, U, h
    zeros(Float64, length(molecules)), 0.0, 0.0) # c_vibr(equilibrium)
end

function compute_mixture!(mixture::Mixture, T::Float64, Tv::Array{Float64,1}, n_atom::Array{Float64,1}, n_molecule::Array{Float64,1};
                          compute_equilibrium=false, use_gupta_yos=false)
    for (i, mol) in enumerate(mixture.molecules)
        mixture.Zr[i] = compute_Z_rot(mol, T)
        mixture.Er[i] = compute_E_rot(mol, T, mixture.Zr[i])

        if compute_equilibrium && !use_gupta_yos
            mixture.Zv[i] = compute_Z_vibr(mol, mixture.vibrational_distribution, T, T)
            mixture.Ev[i] = compute_E_vibr(mol, mixture.vibrational_distribution, T, T, mixture.Zv[i])
            mixture.c_vibr_equilibrium[i] = compute_c_vibr_equilibrium(mol, mixture.vibrational_distribution, T, T, mixture.Zv[i], mixture.Ev[i])
        end
        mixture.Zv[i] = compute_Z_vibr(mol, mixture.vibrational_distribution, T, Tv[i])
        mixture.Ev[i] = compute_E_vibr(mol, mixture.vibrational_distribution, T, Tv[i], mixture.Zv[i])

        mixture.c_rot[i] = compute_c_rot(mol, T, mixture.Zr[i], mixture.Er[i])
        mixture.c_vibrT[i] = compute_c_vibrT(mol, mixture.vibrational_distribution, T, Tv[i], mixture.Zv[i], mixture.Ev[i])
    end

    (mixture.c_v, mixture.c_p) = compute_c_v_and_c_p(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.c_rot, mixture.c_vibrT)
    (mixture.U, mixture.h) = compute_U_and_h(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.Er, mixture.Ev)

    if compute_equilibrium
        if use_gupta_yos
            (mixture.c_v_equilibrium, mixture.c_p_equilibrium) = compute_c_v_and_c_p_equilibrium_gupta_yos(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.c_rot, mixture.c_vibrT)
        else
            (mixture.c_v_equilibrium, mixture.c_p_equilibrium) = compute_c_v_and_c_p(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.c_rot, mixture.c_vibr_equilibrium)
        end
    end
end