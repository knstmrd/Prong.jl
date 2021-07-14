include("particles.jl")
include("constants.jl")

function compute_densities(p, T::Float64, molar_fractions::Array{Float64,1})

    # convert pressure to Pa, 

    return
end

function compute_max_vibr_level(mol::Molecule, T, Tv)
    maxlevel = round(UInt16, (mol.vibrational_energy[2] / (2 * mol.anharmonic_ratio * constants.h * mol.frequency * constants.c)) * (T / Tv) + 0.5)

    if maxlevel > mol.n_vibr
        maxlevel = mol.n_vibr
    end

    return maxlevel
end