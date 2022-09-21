using Unitful
using PhysicalConstants.CODATA2018: k_B
include("utils.jl")


function compute_Z_vibr(mol::Molecule, vd::VibrationalDistribution, T, Tv)
    if mol.anharmonic == false || !vd.use_Treanor
        return sum(exp.(-mol.vibrational_energy ./ (k_B * Tv)))
    else
        if T >= Tv
            return sum(exp.(-mol.Δvibrational_anharmonic ./ (k_B * T)
                       .- mol.vibrational_levels_mult1 ./ (k_B * Tv)))
        else
            maxlevel = compute_max_vibr_level(mol, T, Tv)

            Z = sum(exp.(-mol.Δvibrational_anharmonic[1:maxlevel] ./ (k_B * T)
            .- mol.vibrational_levels_mult1[1:maxlevel] ./ (k_B * Tv)))

            if vd.continue_Treanor_with_Boltzmann && maxlevel + 1 <= mol.n_vibr
                Z += sum(exp.(-mol.vibrational_energy[maxlevel + 1:end] ./ (k_B * Tv)))
            end
            return Z
        end
    end
end

function compute_Z_rot(mol::Molecule, T)
    return sum(mol.rotational_degeneracies .* exp.(-mol.rotational_energy / (k_B * T))) / mol.rotational_symmetry
end

function compute_xi_vibr(mol::Molecule, vd::VibrationalDistribution, T, Tv, Zv::Float64)
    if mol.anharmonic == false || !vd.use_Treanor
        return exp.(-mol.vibrational_energy ./ (k_B * Tv)) / Zv
    else
        if T >= Tv
            return exp.(-mol.Δvibrational_anharmonic ./ (k_B * T)
                     .- mol.vibrational_levels_mult1 ./ (k_B * Tv)) / Zv
        else
            maxlevel = compute_max_vibr_level(mol, T, Tv)

            if vd.continue_Treanor_with_Boltzmann
                res = zeros(mol.n_vibr)

                res[1:maxlevel] .= exp.(-mol.Δvibrational_anharmonic[1:maxlevel] ./ (k_B * T)
                        .- mol.vibrational_levels_mult1[1:maxlevel] ./ (k_B * Tv))
                if maxlevel + 1 <= mol.n_vibr
                    res[maxlevel + 1:end] .= exp.(-mol.vibrational_energy[maxlevel + 1:end] ./ (k_B * Tv))
                end
                return res ./ Zv
            else
                return exp.(-mol.Δvibrational_anharmonic[1:maxlevel] ./ (k_B * T)
                           .- mol.vibrational_levels_mult1[1:maxlevel] ./ (k_B * Tv)) / Zv
            end
        end
    end
end

function compute_E_vibr(mol::Molecule, vd::VibrationalDistribution, T, Tv, Zv::Float64, x_i::Union{Array{Float64,1},Nothing}=nothing)

    if x_i === nothing
        if mol.anharmonic == false || !vd.use_Treanor
            Ev = sum(mol.vibrational_energy .* exp.(-mol.vibrational_energy ./ (k_B * Tv)))
        else
            if T >= Tv
                Ev = sum(mol.vibrational_energy .* exp.(-mol.Δvibrational_anharmonic ./ (k_B * T)
                        .- mol.vibrational_levels_mult1 ./ (k_B * Tv)))
            else
                maxlevel = compute_max_vibr_level(mol, T, Tv)

                if vd.continue_Treanor_with_Boltzmann
                    Ev = sum(mol.vibrational_energy[1:maxlevel] .* exp.(-mol.Δvibrational_anharmonic[1:maxlevel] ./ (k_B * T)
                            .- mol.vibrational_levels_mult1[1:maxlevel] ./ (k_B * Tv)))
                    if maxlevel + 1 <= mol.n_vibr
                        Ev += sum(mol.vibrational_energy[maxlevel + 1:end] .* exp.(-mol.vibrational_energy[maxlevel + 1:end] ./ (k_B * Tv)))
                    end
                else
                    Ev = sum(mol.vibrational_energy[1:maxlevel] .* exp.(-mol.Δvibrational_anharmonic[1:maxlevel] ./ (k_B * T)
                            .- mol.vibrational_levels_mult1[1:maxlevel] ./ (k_B * Tv)))
                end
            end
        end
        return Ev / Zv
    else
        if mol.anharmonic == false || !vd.use_Treanor
            Ev = sum(mol.vibrational_energy .* x_i)
        else
            if T >= Tv
                Ev = sum(mol.vibrational_energy .* x_i)
            else
                maxlevel = compute_max_vibr_level(mol, T, Tv)

                if vd.continue_Treanor_with_Boltzmann
                    Ev = sum(mol.vibrational_energy[1:maxlevel] .* x_i[1:maxlevel])
                    if maxlevel + 1 <= mol.n_vibr
                        Ev += sum(mol.vibrational_energy[maxlevel + 1:end] .* x_i[maxlevel + 1:end])
                    end
                else
                    Ev = sum(mol.vibrational_energy[1:maxlevel] .* x_i[1:maxlevel])
                end
            end
        end
        return Ev
    end
end

function compute_E_rot(mol::Molecule, T, Zr::Float64)
    return sum(mol.rotational_energy .* mol.rotational_degeneracies .* exp.(-mol.rotational_energy / (k_B * T))) / mol.rotational_symmetry / Zr
end

function compute_c_vibrT(mol::Molecule, vd::VibrationalDistribution, T, Tv, Zv, Ev, x_i::Union{Array{Float64,1},Nothing}=nothing)
    if x_i === nothing
        if mol.anharmonic == false || !vd.use_Treanor
            return 0.0u"J * kg^-1 * K^-1"
        else
            if T >= Tv
                exp_distr = exp.(-mol.Δvibrational_anharmonic ./ (k_B * T) .- mol.vibrational_levels_mult1 ./ (k_B * Tv))

                avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* exp_distr)
                avg_i = sum(mol.vibrational_levels .* exp_distr)
                avg_i_ve = sum(mol.vibrational_levels .* mol.vibrational_energy .* exp_distr)

                Ev2 = Ev
            else
                maxlevel = compute_max_vibr_level(mol, T, Tv)

                exp_distr = exp.(-mol.Δvibrational_anharmonic[1:maxlevel] ./ (k_B * T)
                .- mol.vibrational_levels_mult1[1:maxlevel] ./ (k_B * Tv))

                avg_evib_sq = sum(mol.vibrational_energy[1:maxlevel] .* mol.vibrational_energy[1:maxlevel] .* exp_distr)

                avg_i = sum(mol.vibrational_levels[1:maxlevel] .* exp_distr)
                avg_i_ve = sum(mol.vibrational_levels[1:maxlevel] .* mol.vibrational_energy[1:maxlevel] .* exp_distr)
                
                if vd.continue_Treanor_with_Boltzmann && (maxlevel + 1 <= mol.n_vibr)
                    Ev2 = sum(mol.vibrational_energy[1:maxlevel] .* exp_distr) / Zv
                else
                    Ev2 = Ev
                end
            end
        end

        avg_evib_sq /= Zv
        avg_i /= Zv
        avg_i_ve /= Zv


        return (avg_evib_sq - Ev * Ev2 + mol.vibr_energy1 * avg_i * Ev - mol.vibr_energy1 * avg_i_ve) / (k_B * T^2 * mol.mass)
    else
        if mol.anharmonic == false || !vd.use_Treanor
            return 0.0u"J * kg^-1 * K^-1"
        else
            if T >= Tv
                # exp_distr = exp.(-mol.Δvibrational_anharmonic ./ (k_B * T) .- mol.vibrational_levels_mult1 ./ (k_B * Tv))

                avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* x_i)
                avg_i = sum(mol.vibrational_levels .* x_i)
                avg_i_ve = sum(mol.vibrational_levels .* mol.vibrational_energy .* x_i)

                Ev2 = Ev
            else
                maxlevel = compute_max_vibr_level(mol, T, Tv)

                avg_evib_sq = sum(mol.vibrational_energy[1:maxlevel] .* mol.vibrational_energy[1:maxlevel] .* x_i[1:maxlevel])

                avg_i = sum(mol.vibrational_levels[1:maxlevel] .* x_i[1:maxlevel])
                avg_i_ve = sum(mol.vibrational_levels[1:maxlevel] .* mol.vibrational_energy[1:maxlevel] .* x_i[1:maxlevel])
                
                if vd.continue_Treanor_with_Boltzmann && (maxlevel + 1 <= mol.n_vibr)
                    Ev2 = sum(mol.vibrational_energy[1:maxlevel] .* x_i[1:maxlevel]) / Zv
                else
                    Ev2 = Ev
                end
            end
        end

        return (avg_evib_sq - Ev * Ev2 + mol.vibr_energy1 * avg_i * Ev - mol.vibr_energy1 * avg_i_ve) / (k_B * T^2 * mol.mass)
    end
end


function compute_c_vibrTv(mol::Molecule, vd::VibrationalDistribution, T, Tv, Zv, Ev, x_i::Union{Array{Float64,1},Nothing}=nothing)

    if x_i === nothing
        if mol.anharmonic == false || !vd.use_Treanor
            avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* exp.(- mol.vibrational_energy ./ (k_B * Tv)))
            avg_evib_sq /= Zv
            return (avg_evib_sq - Ev^2) / (k_B * Tv^2 * mol.mass)
        else
            if T >= Tv
                exp_distr = exp.(-mol.Δvibrational_anharmonic ./ (k_B * T) .- mol.vibrational_levels_mult1 ./ (k_B * Tv))

                avg_i = sum(mol.vibrational_levels .* exp_distr)
                avg_i_ve = sum(mol.vibrational_levels .* mol.vibrational_energy .* exp_distr)
            else
                maxlevel = compute_max_vibr_level(mol, T, Tv)
                exp_distr = exp.(-mol.Δvibrational_anharmonic[1:maxlevel] ./ (k_B * T) .- mol.vibrational_levels_mult1[1:maxlevel] ./ (k_B * Tv))

                avg_i = sum(mol.vibrational_levels[1:maxlevel] .* exp_distr)
                avg_i_ve = sum(mol.vibrational_levels[1:maxlevel] .* mol.vibrational_energy[1:maxlevel] .* exp_distr)
                
                if vd.continue_Treanor_with_Boltzmann && (maxlevel + 1 <= mol.n_vibr)
                    avg_i += sum(mol.vibrational_energy[maxlevel + 1:end] .* exp.(-mol.vibrational_energy[maxlevel + 1:end] ./ (k_B * Tv))) / mol.vibr_energy1
                    avg_i_ve += sum(mol.vibrational_energy[maxlevel + 1:end] .* mol.vibrational_energy[maxlevel + 1:end] .* exp.(-mol.vibrational_energy[maxlevel + 1:end] ./ (k_B * Tv))) / mol.vibr_energy1
                end
            end
        end

        avg_i /= Zv
        avg_i_ve /= Zv

        return mol.vibr_energy1 * (avg_i_ve - avg_i * Ev) / (k_B * Tv^2 * mol.mass)
    else
        if mol.anharmonic == false || !vd.use_Treanor
            avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* x_i)
            return (avg_evib_sq - Ev^2) / (k_B * Tv^2 * mol.mass)
        else
            if T >= Tv
                avg_i = sum(mol.vibrational_levels .* x_i)
                avg_i_ve = sum(mol.vibrational_levels .* mol.vibrational_energy .* x_i)
            else
                maxlevel = compute_max_vibr_level(mol, T, Tv)

                avg_i = sum(mol.vibrational_levels[1:maxlevel] .* x_i[1:maxlevel])
                avg_i_ve = sum(mol.vibrational_levels[1:maxlevel] .* mol.vibrational_energy[1:maxlevel] .* x_i[1:maxlevel])
                
                if vd.continue_Treanor_with_Boltzmann && (maxlevel + 1 <= mol.n_vibr)
                    avg_i += sum(mol.vibrational_energy[maxlevel + 1:end] .* x_i[maxlevel + 1:end]) / mol.vibr_energy1
                    avg_i_ve += sum(mol.vibrational_energy[maxlevel + 1:end] .* mol.vibrational_energy[maxlevel + 1:end] .* x_i[maxlevel + 1:end]) / mol.vibr_energy1
                end
            end
        end

        return mol.vibr_energy1 * (avg_i_ve - avg_i * Ev) / (k_B * Tv^2 * mol.mass)
    end
end



function compute_c_vibr_equilibrium(mol::Molecule, vd::VibrationalDistribution, T, Tv, Zv, Ev)
    avg_evib_sq = sum(mol.vibrational_energy .* mol.vibrational_energy .* exp.(- mol.vibrational_energy ./ (k_B * Tv)))
    avg_evib_sq /= Zv
    return (avg_evib_sq - Ev^2) / (k_B * T^2 * mol.mass)
end


function compute_c_rot(mol::Molecule, T, Zr, Er)
    avg_erot_sq = sum(mol.rotational_energy .* mol.rotational_energy .* mol.rotational_degeneracies .* exp.(-mol.rotational_energy / (k_B * T))) / mol.rotational_symmetry / Zr
    return (avg_erot_sq - Er^2) / (k_B * T^2 * mol.mass)
end

function compute_c_v_and_c_p_equilibrium_gupta_yos(T, atom_arr::Array{Atom,1}, mol_arr::Array{Molecule,1},
                                       atom_n_arr::Array{typeof(1.0u"m^-3"),1}, mol_n_arr::Array{typeof(1.0u"m^-3"),1},
                                       crot_arr::Array{typeof(1.0u"J * K^-1 * kg^-1"),1},
                                       cvibrT_arr::Array{typeof(1.0u"J * K^-1 * kg^-1"),1})

    T = ustrip(u"K", T)
    if T < 300
        # TODO: fix!
        return 0.0u"J * kg^-1", 0.0u"J * kg^-1"
    elseif T >= 300 && T <= 1000
        coeff_index = 1
    elseif T > 1000 && T <= 6000
        coeff_index = 2
    elseif T > 6000 && T <= 15000
        coeff_index = 3
    elseif T > 15000 && T <= 25000
        coeff_index = 4
    elseif T > 25000
        coeff_index = 5
    end

    rho = 0.0u"kg * m^-3"

    n_tot = sum(atom_n_arr) + sum(mol_n_arr)
    c_p = 0.0u"J * kg^-1"

    for (atom, n) in zip(atom_arr, atom_n_arr)
        rho += atom.mass * n

        c_p += n * atom.gupta_yos_coefficients[1, coeff_index]
        c_p += n * atom.gupta_yos_coefficients[2, coeff_index] * T
        c_p += n * atom.gupta_yos_coefficients[3, coeff_index] * T^2
        c_p += n * atom.gupta_yos_coefficients[4, coeff_index] * T^3
        c_p += n * atom.gupta_yos_coefficients[5, coeff_index] * T^4
    end

    for (mol, n) in zip(mol_arr, mol_n_arr)
        rho += mol.mass * n

        c_p += n * mol.gupta_yos_coefficients[1, coeff_index]
        c_p += n * mol.gupta_yos_coefficients[2, coeff_index] * T
        c_p += n * mol.gupta_yos_coefficients[3, coeff_index] * T^2
        c_p += n * mol.gupta_yos_coefficients[4, coeff_index] * T^3
        c_p += n * mol.gupta_yos_coefficients[5, coeff_index] * T^4
    end

    R_specific = k_B * n_tot / rho
    molar_mass_mixture = rho * constants.N_A / n_tot 

    c_p *= MolarGasConstant / molar_mass_mixture / n_tot  # convert to Joules from calories

    return c_p - R_specific, c_p
end

function compute_c_v_and_c_p(T,
                     atom_arr::Array{Atom,1}, mol_arr::Array{Molecule,1},
                     atom_n_arr::Array{typeof(1.0u"m^-3"),1}, mol_n_arr::Array{typeof(1.0u"m^-3"),1},
                     crot_arr::Array{typeof(1.0u"J * K^-1 * kg^-1"),1},
                     cvibrT_arr::Array{typeof(1.0u"J * K^-1 * kg^-1"),1},
                     n_tot::Union{typeof(1.0u"m^-3"),Nothing}=nothing,
                     rho::Union{typeof(1.0u"kg * m^-3"),Nothing}=nothing)

    if n_tot === nothing
        n_tot = sum(atom_n_arr) + sum(mol_n_arr)
    end

    if rho === nothing
        rho = 0.0u"kg * m^-3"
        for (atom, n) in zip(atom_arr, atom_n_arr)
            rho += atom.mass * n
        end

        for (mol, n) in zip(mol_arr, mol_n_arr)
            rho += mol.mass * n
        end
    end

    c_v = 0.0u"J * m^-3 * K^-1"

    for (mol, n, crot, cvibr) in zip(mol_arr, mol_n_arr, crot_arr, cvibrT_arr)
        c_v += mol.mass * n * (crot + cvibr)
    end

    R_specific = k_B * n_tot / rho

    c_v_out = c_v / rho
    c_v_out += 1.5 * R_specific

    return c_v_out, c_v_out + R_specific
end


function compute_c_v_and_c_p_equilibrium_gupta_yos(T_in, atom_arr::Array{Atom,1}, mol_arr::Array{Molecule,1},
                                                    atom_n_arr::Array{typeof(1.0u"m^-3"),1}, mol_n_arr::Array{typeof(1.0u"m^-3"),1},
                                                    crot_arr::Array{typeof(1.0u"J * K^-1 * kg^-1"),1},
                                                    cvibrT_arr::Array{typeof(1.0u"J * K^-1 * kg^-1"),1})

    T = ustrip(u"K", T_in)
    if T < 300
        # TODO: fix!
        return 0.0u"J * kg^-1", 0.0u"J * kg^-1"
    elseif T >= 300 && T <= 1000
        coeff_index = 1
    elseif T > 1000 && T <= 6000
        coeff_index = 2
    elseif T > 6000 && T <= 15000
        coeff_index = 3
    elseif T > 15000 && T <= 25000
        coeff_index = 4
    elseif T > 25000
        coeff_index = 5
    end

    rho = 0.0u"kg * m^-3"

    n_tot = sum(atom_n_arr) + sum(mol_n_arr)
    c_p = 0.0u"m^-3"

    for (atom, n) in zip(atom_arr, atom_n_arr)
        rho += atom.mass * n

        c_p += n * atom.gupta_yos_coefficients[1, coeff_index]
        c_p += n * atom.gupta_yos_coefficients[2, coeff_index] * T
        c_p += n * atom.gupta_yos_coefficients[3, coeff_index] * T^2
        c_p += n * atom.gupta_yos_coefficients[4, coeff_index] * T^3
        c_p += n * atom.gupta_yos_coefficients[5, coeff_index] * T^4
    end

    for (mol, n) in zip(mol_arr, mol_n_arr)
        rho += mol.mass * n

        c_p += n * mol.gupta_yos_coefficients[1, coeff_index]
        c_p += n * mol.gupta_yos_coefficients[2, coeff_index] * T
        c_p += n * mol.gupta_yos_coefficients[3, coeff_index] * T^2
        c_p += n * mol.gupta_yos_coefficients[4, coeff_index] * T^3
        c_p += n * mol.gupta_yos_coefficients[5, coeff_index] * T^4
    end

    R_specific = k_B * n_tot / rho
    molar_mass_mixture = rho * AvogadroConstant / n_tot 

    c_p_out = c_p * MolarGasConstant / molar_mass_mixture / n_tot

    return c_p_out - R_specific, c_p_out
end


function compute_U_and_h_equilibrium_gupta_yos(T, atom_arr::Array{Atom,1}, mol_arr::Array{Molecule,1},
    atom_n_arr::Array{typeof(1.0u"m^-3"),1}, mol_n_arr::Array{typeof(1.0u"m^-3"),1},
    Erot_arr::Array{typeof(1.0u"J"),1}, Evibr_arr::Array{typeof(1.0u"J"),1})

    T = ustrip(u"K", T)
    if T < 300
        # TODO: fix!
        return 0.0u"J * kg^-1", 0.0u"J * kg^-1"
    elseif T >= 300 && T <= 1000
        coeff_index = 1
    elseif T > 1000 && T <= 6000
        coeff_index = 2
    elseif T > 6000 && T <= 15000
        coeff_index = 3
    elseif T > 15000 && T <= 25000
        coeff_index = 4
    elseif T > 25000
        coeff_index = 5
    end

    rho = 0.0u"kg * m^-3"

    n_tot = sum(atom_n_arr) + sum(mol_n_arr)
    h = 0.0u"m^-3"

    for (atom, n) in zip(atom_arr, atom_n_arr)
        rho += atom.mass * n

        h += n * atom.gupta_yos_coefficients[1, coeff_index]
        h += n * atom.gupta_yos_coefficients[2, coeff_index] * T / 2
        h += n * atom.gupta_yos_coefficients[3, coeff_index] * T^2 / 3
        h += n * atom.gupta_yos_coefficients[4, coeff_index] * T^3 / 4
        h += n * atom.gupta_yos_coefficients[5, coeff_index] * T^4 / 5
        h += n * atom.gupta_yos_coefficients[6, coeff_index] / T
    end

    for (mol, n) in zip(mol_arr, mol_n_arr)
        rho += mol.mass * n

        h += n * mol.gupta_yos_coefficients[1, coeff_index]
        h += n * mol.gupta_yos_coefficients[2, coeff_index] * T / 2
        h += n * mol.gupta_yos_coefficients[3, coeff_index] * T^2 / 3
        h += n * mol.gupta_yos_coefficients[4, coeff_index] * T^3 / 4
        h += n * mol.gupta_yos_coefficients[5, coeff_index] * T^4 / 5
        h += n * mol.gupta_yos_coefficients[6, coeff_index] / T
    end

    molar_mass_mixture = rho * AvogadroConstant / n_tot 

    h_out = h * MolarGasConstant * T * 1.0u"K" / molar_mass_mixture / n_tot  # convert to Joules from calories

    return h_out - n_tot * k_B * T * 1.0u"K" / rho, h_out
end

function compute_U_and_h(T, atom_arr::Array{Atom,1}, mol_arr::Array{Molecule,1},
                        atom_n_arr::Array{typeof(1.0u"m^-3"),1}, mol_n_arr::Array{typeof(1.0u"m^-3"),1},
                        Erot_arr::Array{typeof(1.0u"J"),1}, Evibr_arr::Array{typeof(1.0u"J"),1},
                        n_tot::Union{typeof(1.0u"m^-3"),Nothing}=nothing,
                        rho::Union{typeof(1.0u"kg * m^-3"),Nothing}=nothing)

    if n_tot === nothing
        n_tot = sum(atom_n_arr) + sum(mol_n_arr)
    end

    if rho === nothing
        rho = 0.0u"kg * m^-3"
        for (atom, n) in zip(atom_arr, atom_n_arr)
            rho += atom.mass * n
        end

        for (mol, n) in zip(mol_arr, mol_n_arr)
            rho += mol.mass * n
        end
    end

    U = 0.0u"J * m^-3"

    for (atom, n) in zip(atom_arr, atom_n_arr)
        U += atom.formation_energy * n
    end

    for (mol, n, Erot, Evibr) in zip(mol_arr, mol_n_arr, Erot_arr, Evibr_arr)
        U += (mol.formation_energy + Erot + Evibr) * n
    end

    U_out = U / rho

    U_out += 1.5 * n_tot * k_B * T / rho

    return U_out, U_out + n_tot * k_B * T / rho
end