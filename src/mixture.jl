include("thermodynamic_properties.jl")

mutable struct Mixture
    atoms::Array{Atom,1}
    molecules::Array{Molecule,1}
    vibrational_distribution::VibrationalDistribution

    Zr::Array{Float64,1}
    Er::Array{typeof(1.0u"J"),1}
    Zv::Array{Float64,1}
    Ev::Array{typeof(1.0u"J"),1}

    c_rot::Array{typeof(1.0u"J * K^-1 * kg^-1"),1}
    c_vibrT::Array{typeof(1.0u"J * K^-1 * kg^-1"),1}
    c_vibrTv::Array{typeof(1.0u"J * K^-1 * kg^-1"),1}

    c_v::typeof(1.0u"J * K^-1 * kg^-1")
    c_p::typeof(1.0u"J * K^-1 * kg^-1")

    U::typeof(1.0u"J * kg^-1")
    h::typeof(1.0u"J * kg^-1")

    c_vibr_equilibrium::Array{typeof(1.0u"J * K^-1 * kg^-1"),1}
    c_v_equilibrium::typeof(1.0u"J * K^-1 * kg^-1")
    c_p_equilibrium::typeof(1.0u"J * K^-1 * kg^-1")

    n::typeof(1.0u"m^-3")  # total number density
    rho::typeof(1.0u"kg * m^-3")  # density
end

function create_mixture(atoms, molecules, vd)
    return Mixture(atoms, molecules, vd,
    zeros(Float64, length(molecules)), zeros(typeof(1.0u"J"), length(molecules)),  # Zr/Er
    zeros(Float64, length(molecules)), zeros(typeof(1.0u"J"), length(molecules)),  # Zv/Ev
    zeros(typeof(1.0u"J * K^-1 * kg^-1"), length(molecules)),
    zeros(typeof(1.0u"J * K^-1 * kg^-1"), length(molecules)),
    zeros(typeof(1.0u"J * K^-1 * kg^-1"), length(molecules)), # c_rot, c_vibrT, c_vibrTv
    0.0u"J * K^-1 * kg^-1", # c_v
    0.0u"J * K^-1 * kg^-1", # c_p
    0.0u"J * kg^-1", 0.0u"J * kg^-1", #U, h
    zeros(typeof(1.0u"J * K^-1 * kg^-1"), length(molecules)), # c_vibr(equilibrium)
    0.0u"J * K^-1 * kg^-1", 0.0u"J * K^-1 * kg^-1", # c_v_equilibrium, c_p_equilibrium
    0.0u"m^-3", 0.0u"kg * m^-3") 
end

function compute_mixture!(mixture::Mixture, T, Tv::Array{typeof(1.0u"K"),1}, n_atom::Array{typeof(1.0u"m^-3"),1}, n_molecule::Array{typeof(1.0u"m^-3"),1};
                          compute_equilibrium=false, use_gupta_yos=false)

    mixture.n = sum(n_atom) + sum(n_molecule)

    mixture.rho = 0.0u"kg * m^-3"

    for (i, atom) in enumerate(mixture.atoms)
        mixture.rho += n_atom[i] * atom.mass
    end
    for (i, mol) in enumerate(mixture.molecules)
        mixture.rho += n_molecule[i] * mol.mass
    end

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

    (mixture.c_v, mixture.c_p) = compute_c_v_and_c_p(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.c_rot, mixture.c_vibrT, mixture.n, mixture.rho)
    (mixture.U, mixture.h) = compute_U_and_h(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.Er, mixture.Ev, mixture.n, mixture.rho)

    if compute_equilibrium
        if use_gupta_yos
            (mixture.c_v_equilibrium, mixture.c_p_equilibrium) = compute_c_v_and_c_p_equilibrium_gupta_yos(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.c_rot, mixture.c_vibrT)
        else
            (mixture.c_v_equilibrium, mixture.c_p_equilibrium) = compute_c_v_and_c_p(T, mixture.atoms, mixture.molecules, n_atom, n_molecule, mixture.c_rot, mixture.c_vibr_equilibrium, mixture.n, mixture.rho)
        end
    end
end