include("utils.jl")


function compute_Z_vibr(mol::Molecule, T::Float64, Tv::Float64)
    if mol.anharmonic == false
        return sum(exp.(-mol.vibrational_energy ./ (constants.k * Tv)))
    else
        if T >= Tv
            return sum(exp.(-(mol.vibrational_energy .- mol.vibrational_levels_mult1) ./ (constants.k * Tv)
                       .- mol.vibrational_levels_mult1 ./ (constants.k * T)))
        else
            maxlevel = compute_max_vibr_level(mol, T, Tv)

            if mol.continue_Treanor_with_Boltzmann
                Z = sum(exp.(-(mol.vibrational_energy[1:maxlevel] .- mol.vibrational_levels_mult1[1:maxlevel]) ./ (constants.k * Tv)
                        .- mol.vibrational_levels_mult1[1:maxlevel] ./ (constants.k * T)))
                if maxlevel + 1 <= mol.n_vibr
                    Z += sum(exp.(-(mol.vibrational_energy[maxlevel + 1:end] .- mol.vibrational_levels_mult1[maxlevel + 1:end]) ./ (constants.k * Tv)
                             .- mol.vibrational_levels_mult1[maxlevel + 1:end] ./ (constants.k * T)))
                end
                return Z
            else
                return sum(exp.(-(mol.vibrational_energy[1:maxlevel] .- mol.vibrational_levels_mult1[1:maxlevel]) ./ (constants.k * Tv)
                           .- mol.vibrational_levels_mult1[1:maxlevel] ./ (constants.k * T)))
            end
        end
    end
end

function compute_Z_rot(mol::Molecule, T::Float64)
    return sum(mol.rotational_degeneracies .* exp.(-mol.rotational_energy / (constants.k * T))) / mol.rotational_symmetry
end

function compute_E_vibr(mol::Molecule, T::Float64, Tv::Float64, Zv::Float64)
    if mol.anharmonic == false
        Ev = sum(mol.vibrational_energy .* exp.(-mol.vibrational_energy ./ (constants.k * Tv)))
    else
        if T >= Tv
            Ev = sum(mol.vibrational_energy .* exp.(-(mol.vibrational_energy .- mol.vibrational_levels_mult1) ./ (constants.k * Tv)
                     .- mol.vibrational_levels_mult1 ./ (constants.k * T)))
        else
            maxlevel = compute_max_vibr_level(mol, T, Tv)

            if mol.continue_Treanor_with_Boltzmann
                Ev = sum(mol.vibrational_energy .* exp.(-(mol.vibrational_energy[1:maxlevel] .- mol.vibrational_levels_mult1[1:maxlevel]) ./ (constants.k * Tv)
                        .- mol.vibrational_levels_mult1[1:maxlevel] ./ (constants.k * T)))
                if maxlevel + 1 <= mol.n_vibr
                    Ev += sum(mol.vibrational_energy .* exp.(-(mol.vibrational_energy[maxlevel + 1:end] .- mol.vibrational_levels_mult1[maxlevel + 1:end]) ./ (constants.k * Tv)
                             .- mol.vibrational_levels_mult1[maxlevel + 1:end] ./ (constants.k * T)))
                end
            else
                Ev = sum(mol.vibrational_energy .* exp.(-(mol.vibrational_energy[1:maxlevel] .- mol.vibrational_levels_mult1[1:maxlevel]) ./ (constants.k * Tv)
                           .- mol.vibrational_levels_mult1[1:maxlevel] ./ (constants.k * T)))
            end
        end
    end

    return Ev / Zv
end

function compute_E_rot(mol::Molecule, T::Float64, Zr::Float64)
    return sum(mol.rotational_energy .* mol.rotational_degeneracies .* exp.(-mol.rotational_energy / (constants.k * T))) / mol.rotational_symmetry / Zr
end

function compute_c_vibr(mol::Molecule, T::Float64, Tv::Float64, Zv::Float64, Ev::Float64)
    if mol.anharmonic == false
        avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* exp.(-mol.vibrational_energy ./ (constants.k * Tv)))
    else
        if T >= Tv
            avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* exp.(-(mol.vibrational_energy .- mol.vibrational_levels_mult1) ./ (constants.k * Tv)
                     .- mol.vibrational_levels_mult1 ./ (constants.k * T)))
        else
            maxlevel = compute_max_vibr_level(mol, T, Tv)

            if mol.continue_Treanor_with_Boltzmann
                avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* exp.(-(mol.vibrational_energy[1:maxlevel] .- mol.vibrational_levels_mult1[1:maxlevel]) ./ (constants.k * Tv)
                        .- mol.vibrational_levels_mult1[1:maxlevel] ./ (constants.k * T)))
                if maxlevel + 1 <= mol.n_vibr
                    avg_evib_sq += sum(mol.vibrational_energy .* mol.vibrational_energy .* exp.(-(mol.vibrational_energy[maxlevel + 1:end] .- mol.vibrational_levels_mult1[maxlevel + 1:end]) ./ (constants.k * Tv)
                             .- mol.vibrational_levels_mult1[maxlevel + 1:end] ./ (constants.k * T)))
                end
            else
                avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* exp.(-(mol.vibrational_energy[1:maxlevel] .- mol.vibrational_levels_mult1[1:maxlevel]) ./ (constants.k * Tv)
                           .- mol.vibrational_levels_mult1[1:maxlevel] ./ (constants.k * T)))
            end
        end
    end

    avg_evib_sq /= Zv
    return (avg_evib_sq - Ev^2) / (constants.k * T^2 * mol.mass)
end

function compute_c_rot(mol::Molecule, T::Float64, Zr::Float64, Er::Float64)
    avg_erot_sq = sum(mol.rotational_energy .* mol.rotational_energy .* mol.rotational_degeneracies .* exp.(-mol.rotational_energy / (constants.k * T))) / mol.rotational_symmetry / Zr
    return (avg_erot_sq - Er^2) / (constants.k * T^2 * mol.mass)
end

function compute_c_v_and_c_p(T::Float64,
                     atom_arr::Array{Molecule,1}, mol_arr::Array{Molecule,1},
                     atom_n_arr::Array{Float64,1}, mol_n_arr::Array{Float64,1},
                     crot_arr::Array{Float64,1}, cvibr_arr::Array{Float64,1})

        rho = 0.0

        n_tot = sum(atom_n_arr) + sum(mol_n_arr)
        c_v = 0.0

        for (atom, n) in zip(atom_arr, atom_n_arr)
            rho += atom.mass * n
        end

        for (mol, n, crot, cvibr) in zip(mol_arr, atom_n_arr, crot_arr, cvibr_arr)
            rho += mol.mass * n
            c_v += mol.mass * n * (crot + cvibr)
        end

        R_specific = constants.k * n_tot / rho

        c_v /= rho
        c_v += 1.5 * R_specific

        return c_v, c_v + R_specific
end


function compute_U_and_h(T::Float64, atom_arr::Array{Molecule,1}, mol_arr::Array{Molecule,1},
                   atom_n_arr::Array{Float64,1}, mol_n_arr::Array{Float64,1},
                   Erot_arr::Array{Float64,1}, Evibr_arr::Array{Float64,1})
    U = 0.0

    rho = 0.0

    n_tot = sum(atom_n_arr) + sum(mol_n_arr)

    for (atom, n) in zip(atom_arr, atom_n_arr)
        rho += atom.mass * n
        U += atom.formation_energy * n
    end

    for (mol, n, Erot, Evibr) in zip(mol_arr, atom_n_arr, Erot_arr, Evibr_arr)
        rho += mol.mass * n
        U += (mol.formation_energy + Erot + Evibr) * n
    end

    U /= rho

    U += 1.5 * n_tot * constants.k * T / rho

    return c_v, c_v + R_specific
end