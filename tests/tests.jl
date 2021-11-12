using YAML: length
include("../src/mixture.jl")
using Test

global const molnames = ["N2", "O2", "NO"]

@testset "data loading" begin
    rtol = 0.0001

    for molname in molnames
        mol = create_molecule("data/particles.yaml", molname, simplified_anharmonic=true)
        @test last(mol.vibrational_energy) < mol.dissociation_energy
        @test mol.n_vibr - 1 == UInt16(mol.vibrational_levels[length(mol.vibrational_levels)])

        molh = create_molecule("data/particles.yaml", molname, anharmonic=false)
        @test last(molh.vibrational_energy) < molh.dissociation_energy
        @test molh.n_vibr - 1 == UInt16(molh.vibrational_levels[length(molh.vibrational_levels)])

        @test molh.n_vibr < mol.n_vibr

        @test true == isapprox(molh.vibrational_energy[2], molh.vibrational_levels_mult1[2], rtol=rtol)
        @test true == isapprox(molh.vibrational_energy[molh.n_vibr], molh.vibrational_levels_mult1[molh.n_vibr], rtol=rtol)
    end
end


@testset "distributions" begin
    rtol = 0.0001

    for molname in molnames
        mol = create_molecule("data/particles.yaml", molname, anharmonic=false)

        vd = create_vibrational_distribution(false, false, false)
        for T in [500.0u"K", 5000.0u"K", 20000.0u"K"]
            for Tv in [500.0u"K", 5000.0u"K", 20000.0u"K"]
                Zv = compute_Z_vibr(mol, vd, T, Tv)

                xi_arr = compute_xi_vibr(mol, vd, T, Tv, Zv)

                @test length(xi_arr) == mol.n_vibr
                @test true == isapprox(sum(xi_arr), 1.0, rtol=rtol)
                @test true == all(xi_arr[1:length(xi_arr) - 1] .>= xi_arr[2:length(xi_arr)])
            end
        end
    end

    for molname in molnames
        mol = create_molecule("data/particles.yaml", molname, anharmonic=true, simplified_anharmonic=true)
        vd_boltzmann = create_vibrational_distribution(true, false, false)
        vd_tr = create_vibrational_distribution(true, true, false)
        vd_tr_cont = create_vibrational_distribution(true, true, true)

        for vd in [vd_boltzmann, vd_tr, vd_tr_cont]
            for T in [500.0u"K", 5000.0u"K", 20000.0u"K"]
                for Tv in [500.0u"K", 5000.0u"K", 20000.0u"K"]
                    Zv = compute_Z_vibr(mol, vd, T, Tv)

                    xi_arr = compute_xi_vibr(mol, vd, T, Tv, Zv)

                    if vd != vd_tr
                        @test length(xi_arr) == mol.n_vibr
                    else
                        @test length(xi_arr) <= mol.n_vibr
                    end
                    @test true == isapprox(sum(xi_arr), 1.0, rtol=rtol)
                    @test true == all(xi_arr[1:length(xi_arr) - 1] .>= xi_arr[2:length(xi_arr)])
                end
            end
        end
    end
end

@testset "specific heats" begin    

    # check that c_rot is approximately k/m, ±1%

    rtol = 0.01
    ΔT = 1u"K"

    for molname in molnames
        mol = create_molecule("data/particles.yaml", molname, anharmonic=true, simplified_anharmonic=true)
        for T in [500.0u"K", 1000.0u"K", 2500.0u"K", 5000.0u"K"]
            Zr = compute_Z_rot(mol, T)
            Er = compute_E_rot(mol, T, Zr)
            crot = compute_c_rot(mol, T, Zr, Er)
            crot_const = k_B / mol.mass
            @test true == isapprox(crot, crot_const, rtol=rtol)


            Zrm = compute_Z_rot(mol, T - ΔT)
            Erm = compute_E_rot(mol, T - ΔT, Zrm) / mol.mass
            Zrp = compute_Z_rot(mol, T + ΔT)
            Erp = compute_E_rot(mol, T + ΔT, Zrp) / mol.mass

            @test true == isapprox(crot, (Erp - Erm) / (2 * ΔT), rtol=rtol)
        end

        # crot is not very close to k/m at higher temperatures, so we're dropping that test
        for T in [7500.0u"K", 10000.0u"K", 15000.0u"K", 20000.0u"K"]
            Zr = compute_Z_rot(mol, T)
            Er = compute_E_rot(mol, T, Zr)
            crot = compute_c_rot(mol, T, Zr, Er)

            Zrm = compute_Z_rot(mol, T - ΔT)
            Erm = compute_E_rot(mol, T - ΔT, Zrm) / mol.mass
            Zrp = compute_Z_rot(mol, T + ΔT)
            Erp = compute_E_rot(mol, T + ΔT, Zrp) / mol.mass

            @test true == isapprox(crot, (Erp - Erm) / (2 * ΔT), rtol=rtol)
        end
    end

    for molname in molnames
        mol = create_molecule("data/particles.yaml", molname, anharmonic=false)

        vd = create_vibrational_distribution(false, false, false)

        for T in [500.0u"K", 5000.0u"K", 20000.0u"K"]
            for Tv in [500.0u"K", 5000.0u"K", 20000.0u"K"]
                Zv = compute_Z_vibr(mol, vd, T, Tv)
                Ev = compute_E_vibr(mol, vd, T, Tv, Zv)
                c_vibrT = compute_c_vibrT(mol, vd, T, Tv, Zv, Ev)
                c_vibrTv = compute_c_vibrTv(mol, vd, T, Tv, Zv, Ev)


                Zvm = compute_Z_vibr(mol, vd, T - ΔT, Tv)
                Evm = compute_E_vibr(mol, vd, T - ΔT, Tv, Zvm) / mol.mass
                Zvp = compute_Z_vibr(mol, vd, T + ΔT, Tv)
                Evp = compute_E_vibr(mol, vd, T + ΔT, Tv, Zvp) / mol.mass

                @test c_vibrT <= 0.0u"J * kg^-1 * K^-1"
                @test true == isapprox(c_vibrT, (Evp - Evm) / (2 * ΔT), rtol=rtol)


                Zvm = compute_Z_vibr(mol, vd, T, Tv - ΔT)
                Evm = compute_E_vibr(mol, vd, T, Tv - ΔT, Zvm) / mol.mass
                Zvp = compute_Z_vibr(mol, vd, T, Tv + ΔT)
                Evp = compute_E_vibr(mol, vd, T, Tv + ΔT, Zvp) / mol.mass

                @test true == isapprox(c_vibrTv, (Evp - Evm) / (2 * ΔT), rtol=rtol)
            end
        end
    end


    # ΔT = 0.05
    for molname in molnames
        mol = create_molecule("data/particles.yaml", molname, anharmonic=true, simplified_anharmonic=true)

        for continue_Treanor_with_Boltzmann in [true, false]
            # continue_Treanor_with_Boltzmann=continue_Treanor_with_Boltzmann)
            vd = create_vibrational_distribution(true, true, continue_Treanor_with_Boltzmann)
            for T in [500.0u"K", 5000.0u"K", 20000.0u"K"]
                for Tv in [500.0u"K", 5000.0u"K", 20000.0u"K"]
                    Zv = compute_Z_vibr(mol, vd, T, Tv)
                    Ev = compute_E_vibr(mol, vd, T, Tv, Zv)
                    c_vibrT = compute_c_vibrT(mol, vd, T, Tv, Zv, Ev)
                    c_vibrTv = compute_c_vibrTv(mol, vd, T, Tv, Zv, Ev)


                    Zvm = compute_Z_vibr(mol, vd, T - ΔT, Tv)
                    Evm = compute_E_vibr(mol, vd, T - ΔT, Tv, Zvm) / mol.mass
                    Zvp = compute_Z_vibr(mol, vd, T + ΔT, Tv)
                    Evp = compute_E_vibr(mol, vd, T + ΔT, Tv, Zvp) / mol.mass

                    # some numerical precision issues arise for this distribution, so we relax the tolerances a bit
                    if (abs(c_vibrT) > 1e-8u"J * kg^-1 * K^-1")
                        @test true == isapprox(c_vibrT, (Evp - Evm) / (2 * ΔT), rtol=rtol)
                    else
                        @test true == isapprox(c_vibrT, (Evp - Evm) / (2 * ΔT), atol=1e-10u"J * kg^-1 * K^-1")
                    end

                    Zvm = compute_Z_vibr(mol, vd, T, Tv - ΔT)
                    Evm = compute_E_vibr(mol, vd, T, Tv - ΔT, Zvm) / mol.mass
                    Zvp = compute_Z_vibr(mol, vd, T, Tv + ΔT)
                    Evp = compute_E_vibr(mol, vd, T, Tv + ΔT, Zvp) / mol.mass

                    @test true == isapprox(c_vibrTv, (Evp - Evm) / (2 * ΔT), rtol=rtol)
                end
            end
        end

        vd = create_vibrational_distribution(true, false, false)
        for T in [500.0u"K", 5000.0u"K", 20000.0u"K"]
            for Tv in [500.0u"K", 5000.0u"K", 20000.0u"K"]
                Zv = compute_Z_vibr(mol, vd, T, Tv)
                Ev = compute_E_vibr(mol, vd, T, Tv, Zv)
                c_vibrT = compute_c_vibrT(mol, vd, T, Tv, Zv, Ev)
                c_vibrTv = compute_c_vibrTv(mol, vd, T, Tv, Zv, Ev)


                Zvm = compute_Z_vibr(mol, vd, T - ΔT, Tv)
                Evm = compute_E_vibr(mol, vd, T - ΔT, Tv, Zvm) / mol.mass
                Zvp = compute_Z_vibr(mol, vd, T + ΔT, Tv)
                Evp = compute_E_vibr(mol, vd, T + ΔT, Tv, Zvp) / mol.mass

                @test c_vibrT <= 0.0u"J * kg^-1 * K^-1"
                @test true == isapprox(c_vibrT, (Evp - Evm) / (2 * ΔT), rtol=rtol)

                Zvm = compute_Z_vibr(mol, vd, T, Tv - ΔT)
                Evm = compute_E_vibr(mol, vd, T, Tv - ΔT, Zvm) / mol.mass
                Zvp = compute_Z_vibr(mol, vd, T, Tv + ΔT)
                Evp = compute_E_vibr(mol, vd, T, Tv + ΔT, Zvp) / mol.mass

                @test true == isapprox(c_vibrTv, (Evp - Evm) / (2 * ΔT), rtol=rtol)
            end
        end
    end
end


@testset "specific heats curve fits" begin

    T_arr = [500.0u"K", 1000.0u"K", 2500.0u"K"]
    rtol = 0.03

    for atomname in ["N", "O"]
        for T in T_arr
            atom = create_atom("data/particles.yaml", atomname)

            molarr = Molecule[]
            nmolarr = typeof(1.0u"m^-3")[]
            crotarr = typeof(1.0u"J * kg^-1 * K^-1")[]
            cvibrarr = typeof(1.0u"J * kg^-1 * K^-1")[]

            c_v_c_p_gy = compute_c_v_and_c_p_equilibrium_gupta_yos(T, [atom], molarr, [1e25u"m^-3"], nmolarr, crotarr, cvibrarr)
            c_v_c_p = compute_c_v_and_c_p(T, [atom], molarr, [1e25u"m^-3"], nmolarr, crotarr, cvibrarr)
            @test true == isapprox(c_v_c_p[2], c_v_c_p_gy[2], rtol=rtol)
        end
    end


    T_arr = [500.0u"K", 1000.0u"K", 2000.0u"K"]
    for molname in molnames
        mol = create_molecule("data/particles.yaml", molname, anharmonic=true, simplified_anharmonic=true)
        vd = create_vibrational_distribution(true, true, false)
        for T in T_arr
            atomarr = Atom[]
            natomarr = typeof(1.0u"m^-3")[]
            crotarr = typeof(1.0u"J * kg^-1 * K^-1")[]
            cvibrarr = typeof(1.0u"J * kg^-1 * K^-1")[]

            c_v_c_p_gy = compute_c_v_and_c_p_equilibrium_gupta_yos(T, atomarr, [mol], natomarr, [1e23u"m^-3"], crotarr, cvibrarr)

            Zr = compute_Z_rot(mol, T)
            Er = compute_E_rot(mol, T, Zr)
            crotarr = [compute_c_rot(mol, T, Zr, Er)]
            Zv = compute_Z_vibr(mol, vd, T, T)
            Ev = compute_E_vibr(mol, vd, T, T, Zv)
            crotarr = [compute_c_rot(mol, T, Zr, Er)]
            cvibrarr = [compute_c_vibrTv(mol, vd, T, T, Zv, Ev)]

            c_v_c_p = compute_c_v_and_c_p(T, atomarr, [mol], natomarr, [1e23u"m^-3"], crotarr, cvibrarr)
            @test true == isapprox(c_v_c_p[2], c_v_c_p_gy[2], rtol=rtol)
            # println("T: ", T, ", ", c_v_c_p[1])
        end
    end
end

@testset "enthalpy curve fits" begin
    T_arr = [500.0u"K", 1500.0u"K", 2000.0u"K"]
    rtol = 0.04

    ΔT = 300u"K"
    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname, anharmonic=true)
        vd = create_vibrational_distribution(true, true, false)
        for T in T_arr

            atomarr = Atom[]
            natomarr = typeof(1.0u"m^-3")[]
            crotarr = typeof(1.0u"J * kg^-1 * K^-1")[]
            cvibrarr = typeof(1.0u"J * kg^-1 * K^-1")[]

            Zr = compute_Z_rot(mol, T)
            Er_arr = [compute_E_rot(mol, T, Zr)]
            Zv = compute_Z_vibr(mol, vd, T, T)
            Ev_arr = [compute_E_vibr(mol, vd, T, T, Zv)]

            U_h_1 = compute_U_and_h(T, atomarr, [mol], natomarr, [1e23u"m^-3"], Er_arr, Ev_arr)
            U_h_1_gy = compute_U_and_h_equilibrium_gupta_yos(T, atomarr, [mol], natomarr, [1e23u"m^-3"], Er_arr, Ev_arr)

            Zr = compute_Z_rot(mol, T + ΔT)
            Er_arr = [compute_E_rot(mol, T + ΔT, Zr)]
            Zv = compute_Z_vibr(mol, vd, T + ΔT, T + ΔT)
            Ev_arr = [compute_E_vibr(mol, vd, T + ΔT, T + ΔT, Zv)]

            U_h_2 = compute_U_and_h(T + ΔT, atomarr, [mol], natomarr, [1e23u"m^-3"], Er_arr, Ev_arr)
            U_h_2_gy = compute_U_and_h_equilibrium_gupta_yos(T + ΔT, atomarr, [mol], natomarr, [1e23u"m^-3"], Er_arr, Ev_arr)


            @test true == isapprox(U_h_2[2] - U_h_1[2], U_h_2_gy[2] - U_h_1_gy[2], rtol=rtol)
        end
    end
end

@testset "mixture specific heats" begin

    T_arr = [500.0u"K", 5000.0u"K", 10000.0u"K"]
    rtol = 0.01
    ΔT = 0.05u"K"
    n = 1e23u"m^-3"

    rtol_mix = 0.00001

    for (atomname, molname) in zip(["N", "O", "N"], ["N2", "O2", "NO"])
        atom = create_atom("data/particles.yaml", atomname)

        # in a Treanor distribution it's more complicated, because both T and Tv appear in the distribution
        # for a harmonic oscillator everything's nice though, and we're testing thermal equilibrium here
        mol = create_molecule("data/particles.yaml", molname, anharmonic=false)

        vd = create_vibrational_distribution(false, false, false)

        mixture = create_mixture([atom], [mol], vd)
        for T in T_arr

            Zrm = compute_Z_rot(mol, T - ΔT)
            Er_arrm = [compute_E_rot(mol, T - ΔT, Zrm)]
            Zvm = compute_Z_vibr(mol, vd, T - ΔT, T - ΔT)
            Ev_arrm = [compute_E_vibr(mol, vd, T - ΔT, T - ΔT, Zvm)]

            Zrp = compute_Z_rot(mol, T + ΔT)
            Er_arrp = [compute_E_rot(mol, T + ΔT, Zrp)]
            Zvp = compute_Z_vibr(mol, vd, T + ΔT, T + ΔT)
            Ev_arrp = [compute_E_vibr(mol, vd, T + ΔT, T + ΔT, Zvp)]

            Zr = compute_Z_rot(mol, T)
            Er = compute_E_rot(mol, T, Zr)
            Zv = compute_Z_vibr(mol, vd, T, T)
            Ev = compute_E_vibr(mol, vd, T, T, Zv)

            crotarr = [compute_c_rot(mol, T, Zr, Er)]
            cvibrarr = [compute_c_vibrTv(mol, vd, T, T, Zv, Ev)]

            for x_atom in [0.0, 0.25, 0.5, 1.0]
                U_h_1 = compute_U_and_h(T - ΔT, [atom], [mol], [n * x_atom], [(1.0 - x_atom) * n], Er_arrm, Ev_arrm)
                U_h_2 = compute_U_and_h(T + ΔT, [atom], [mol], [n * x_atom], [(1.0 - x_atom) * n], Er_arrp, Ev_arrp)
                
                c_v_c_p = compute_c_v_and_c_p(T, [atom], [mol], [n * x_atom], [(1.0 - x_atom) * n], crotarr, cvibrarr)

                compute_mixture!(mixture, T - ΔT, [T - ΔT], [n * x_atom], [(1.0 - x_atom) * n])
                @test true == isapprox(c_v_c_p[1], (U_h_2[1] - U_h_1[1]) / (2 * ΔT), rtol=rtol)
                @test true == isapprox(mixture.U, U_h_1[1], rtol=rtol_mix)
            end
        end
    end
end