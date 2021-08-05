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
end

function create_mixture(atoms, molecules, vd)
    return Mixture(atoms, molecules, vd,
    zeros(Float64, length(molecules)), zeros(Float64, length(molecules)),  # Zr/Er
    zeros(Float64, length(molecules)), zeros(Float64, length(molecules)),  # Zv/Ev
    zeros(Float64, length(molecules)), zeros(Float64, length(molecules)), zeros(Float64, length(molecules)), # c_rot, c_vibrT, c_vibrTv
    0.0, 0.0, 0.0, 0.0)  # c_v, c_p, U, h
end

function compute_mixture!(mixture::Mixture, T::Float64, Tv::Array{Float64,1}, n_atom::Array{Float64,1}, n_molecule::Array{Float64,1})
    for (i, mol) in enumerate(mixture.molecules)
        mixture.Zr[i] = compute_Z_rot(mol, T)
        mixture.Er[i] = compute_E_rot(mol, T, mixture.Zr[i])

        mixture.Zv[i] = compute_Z_vibr(mol, mixture.vibrational_distribution, T, Tv[i])
        mixture.Ev[i] = compute_E_vibr(mol, mixture.vibrational_distribution, T, Tv[i], mixture.Zv[i])

        mixture.c_rot[i] = compute_c_rot(mol, T, mixture.Zr[i], mixture.Er[i])
        mixture.c_vibrT[i] = compute_c_vibrT(mol, mixture.vibrational_distribution, T, Tv[i], mixture.Zv[i], mixture.Ev[i])

        (mixture.c_v, mixture.c_p) = compute_c_v_and_c_p(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.c_rot, mixture.c_vibrT)
        (mixture.U, mixture.h) = compute_U_and_h(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.Er, mixture.Ev)
    end
end