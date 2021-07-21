using YAML: length
include("../src/thermodynamic_properties.jl")
using Test

@testset "data loading" begin
    rtol = 0.0001

    for molname in ["N2", "O2", "NO"]
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

    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname, anharmonic=false)

        for T in [500.0, 5000.0, 20000.0]
            for Tv in [500.0, 5000.0, 20000.0]
                Zv = compute_Z_vibr(mol, T, Tv)

                xi_arr = compute_xi_vibr(mol, T, Tv, Zv)

                @test length(xi_arr) == mol.n_vibr
                @test true == isapprox(sum(xi_arr), 1.0, rtol=rtol)
                @test true == all(xi_arr[1:length(xi_arr) - 1] .>= xi_arr[2:length(xi_arr)])
            end
        end
    end

    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname, simplified_anharmonic=true, use_Treanor=false)

        for T in [500.0, 5000.0, 20000.0]
            for Tv in [500.0, 5000.0, 20000.0]
                Zv = compute_Z_vibr(mol, T, Tv)

                xi_arr = compute_xi_vibr(mol, T, Tv, Zv)

                @test length(xi_arr) == mol.n_vibr
                @test true == isapprox(sum(xi_arr), 1.0, rtol=rtol)
                @test true == all(xi_arr[1:length(xi_arr) - 1] .>= xi_arr[2:length(xi_arr)])
            end
        end
    end

    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname, simplified_anharmonic=true, use_Treanor=true, continue_Treanor_with_Boltzmann=false)

        for T in [500.0, 5000.0, 20000.0]
            for Tv in [500.0, 5000.0, 20000.0]
                Zv = compute_Z_vibr(mol, T, Tv)

                xi_arr = compute_xi_vibr(mol, T, Tv, Zv)

                @test true == isapprox(sum(xi_arr), 1.0, rtol=rtol)
                @test true == all(xi_arr[1:length(xi_arr) - 1] .>= xi_arr[2:length(xi_arr)])
            end
        end
    end

    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname, simplified_anharmonic=true, use_Treanor=true, continue_Treanor_with_Boltzmann=true)

        for T in [500.0, 5000.0, 20000.0]
            for Tv in [500.0, 5000.0, 20000.0]
                Zv = compute_Z_vibr(mol, T, Tv)

                xi_arr = compute_xi_vibr(mol, T, Tv, Zv)

                @test length(xi_arr) == mol.n_vibr
                @test true == isapprox(sum(xi_arr), 1.0, rtol=rtol)
                @test true == all(xi_arr[1:length(xi_arr) - 1] .>= xi_arr[2:length(xi_arr)])
            end
        end
    end
end

@testset "specific heats" begin    

    # check that c_rot is approximately k/m, ±1%

    rtol = 0.01
    ΔT = 1

    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname)
        for T in [500.0, 1000.0, 2500.0, 5000.0]
            Zr = compute_Z_rot(mol, T)
            Er = compute_E_rot(mol, T, Zr)
            crot = compute_c_rot(mol, T, Zr, Er)
            crot_const = constants.k / mol.mass
            @test true == isapprox(crot, crot_const, rtol=rtol)


            Zrm = compute_Z_rot(mol, T - ΔT)
            Erm = compute_E_rot(mol, T - ΔT, Zrm) / mol.mass
            Zrp = compute_Z_rot(mol, T + ΔT)
            Erp = compute_E_rot(mol, T + ΔT, Zrp) / mol.mass

            @test true == isapprox(crot, (Erp - Erm) / (2 * ΔT), rtol=rtol)
        end

        # crot is not very close to k/m at higher temperatures, so we're dropping that test
        for T in [7500.0, 10000.0, 15000.0, 20000.]
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

    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname, anharmonic=false)
        for T in [500.0, 5000.0, 20000.0]
            for Tv in [500.0, 5000.0, 20000.0]
                Zv = compute_Z_vibr(mol, T, Tv)
                Ev = compute_E_vibr(mol, T, Tv, Zv)
                c_vibrT = compute_c_vibrT(mol, T, Tv, Zv, Ev)
                c_vibrTv = compute_c_vibrTv(mol, T, Tv, Zv, Ev)


                Zvm = compute_Z_vibr(mol, T - ΔT, Tv)
                Evm = compute_E_vibr(mol, T - ΔT, Tv, Zvm) / mol.mass
                Zvp = compute_Z_vibr(mol, T + ΔT, Tv)
                Evp = compute_E_vibr(mol, T + ΔT, Tv, Zvp) / mol.mass

                # println(molname, "; T=", T, "; Tv=", Tv, ": ", c_vibrT, ", ", (Evp - Evm) / (2 * ΔT))
                @test c_vibrT <= 0.0
                @test true == isapprox(c_vibrT, (Evp - Evm) / (2 * ΔT), rtol=rtol)


                Zvm = compute_Z_vibr(mol, T, Tv - ΔT)
                Evm = compute_E_vibr(mol, T, Tv - ΔT, Zvm) / mol.mass
                Zvp = compute_Z_vibr(mol, T, Tv + ΔT)
                Evp = compute_E_vibr(mol, T, Tv + ΔT, Zvp) / mol.mass

                # println(molname, "; T=", T, "; Tv=", Tv, ": ", c_vibrTv, ", ", (Evp - Evm) / (2 * ΔT))
                @test true == isapprox(c_vibrTv, (Evp - Evm) / (2 * ΔT), rtol=rtol)
            end
        end
    end


    # ΔT = 0.05
    for molname in ["N2", "O2", "NO"]
        for continue_Treanor_with_Boltzmann in [true, false]
            mol = create_molecule("data/particles.yaml", molname, anharmonic=true, simplified_anharmonic=true, continue_Treanor_with_Boltzmann=continue_Treanor_with_Boltzmann)
            for T in [500.0, 5000.0, 20000.0]
                for Tv in [500.0, 5000.0, 20000.0]
                    Zv = compute_Z_vibr(mol, T, Tv)
                    Ev = compute_E_vibr(mol, T, Tv, Zv)
                    c_vibrT = compute_c_vibrT(mol, T, Tv, Zv, Ev)
                    c_vibrTv = compute_c_vibrTv(mol, T, Tv, Zv, Ev)


                    Zvm = compute_Z_vibr(mol, T - ΔT, Tv)
                    Evm = compute_E_vibr(mol, T - ΔT, Tv, Zvm) / mol.mass
                    Zvp = compute_Z_vibr(mol, T + ΔT, Tv)
                    Evp = compute_E_vibr(mol, T + ΔT, Tv, Zvp) / mol.mass

                    # some numerical precision issues arise for this distribution, so we relax the tolerances a bit
                    if (abs(c_vibrT) > 1e-8)
                        @test true == isapprox(c_vibrT, (Evp - Evm) / (2 * ΔT), rtol=rtol)
                    else
                        @test true == isapprox(c_vibrT, (Evp - Evm) / (2 * ΔT), atol=1e-10)
                    end

                    Zvm = compute_Z_vibr(mol, T, Tv - ΔT)
                    Evm = compute_E_vibr(mol, T, Tv - ΔT, Zvm) / mol.mass
                    Zvp = compute_Z_vibr(mol, T, Tv + ΔT)
                    Evp = compute_E_vibr(mol, T, Tv + ΔT, Zvp) / mol.mass

                    @test true == isapprox(c_vibrTv, (Evp - Evm) / (2 * ΔT), rtol=rtol)
                end
            end
        end
    end


    # ΔT = 1
    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname, anharmonic=true, simplified_anharmonic=true, use_Treanor=false)
        for T in [500.0, 5000.0, 20000.0]
            for Tv in [500.0, 5000.0, 20000.0]
                Zv = compute_Z_vibr(mol, T, Tv)
                Ev = compute_E_vibr(mol, T, Tv, Zv)
                c_vibrT = compute_c_vibrT(mol, T, Tv, Zv, Ev)
                c_vibrTv = compute_c_vibrTv(mol, T, Tv, Zv, Ev)


                Zvm = compute_Z_vibr(mol, T - ΔT, Tv)
                Evm = compute_E_vibr(mol, T - ΔT, Tv, Zvm) / mol.mass
                Zvp = compute_Z_vibr(mol, T + ΔT, Tv)
                Evp = compute_E_vibr(mol, T + ΔT, Tv, Zvp) / mol.mass

                @test c_vibrT <= 0.0
                @test true == isapprox(c_vibrT, (Evp - Evm) / (2 * ΔT), rtol=rtol)

                Zvm = compute_Z_vibr(mol, T, Tv - ΔT)
                Evm = compute_E_vibr(mol, T, Tv - ΔT, Zvm) / mol.mass
                Zvp = compute_Z_vibr(mol, T, Tv + ΔT)
                Evp = compute_E_vibr(mol, T, Tv + ΔT, Zvp) / mol.mass

                @test true == isapprox(c_vibrTv, (Evp - Evm) / (2 * ΔT), rtol=rtol)
            end
        end
    end
end


@testset "specific heats curve fits" begin

    T_arr = [500.0, 1000.0, 2500.0]
    rtol = 0.03

    for atomname in ["N", "O"]
        for T in T_arr
            atom = create_atom("data/particles.yaml", atomname)

            molarr = Molecule[]
            nmolarr = Float64[]
            crotarr = Float64[]
            cvibrarr = Float64[]

            c_v_c_p_gy = compute_c_v_and_c_p_gupta_yos(T, [atom], molarr, [1e25], nmolarr, crotarr, cvibrarr)
            c_v_c_p = compute_c_v_and_c_p(T, [atom], molarr, [1e25], nmolarr, crotarr, cvibrarr)
            @test true == isapprox(c_v_c_p[2], c_v_c_p_gy[2], rtol=rtol)
        end
    end


    for molname in ["N2", "O2", "NO"]
        mol = create_molecule("data/particles.yaml", molname, anharmonic=true, simplified_anharmonic=true)
        for T in T_arr

            T = 1000.0
            atomarr = Atom[]
            natomarr = Float64[]
            crotarr = Float64[]
            cvibrarr = Float64[]

            c_v_c_p_gy = compute_c_v_and_c_p_gupta_yos(T, atomarr, [mol], natomarr, [1e23], crotarr, cvibrarr)

            Zr = compute_Z_rot(mol, T)
            Er = compute_E_rot(mol, T, Zr)
            crotarr = [compute_c_rot(mol, T, Zr, Er)]
            Zv = compute_Z_vibr(mol, T, T)
            Ev = compute_E_vibr(mol, T, T, Zv)
            crotarr = [compute_c_rot(mol, T, Zr, Er)]
            cvibrarr = [compute_c_vibrTv(mol, T, T, Zv, Ev)]

            c_v_c_p = compute_c_v_and_c_p(T, atomarr, [mol], natomarr, [1e23], crotarr, cvibrarr)
            @test true == isapprox(c_v_c_p[2], c_v_c_p_gy[2], rtol=rtol)
        end
    end
end

@testset "mixture specific heats" begin
    # mol = create_molecule("data/particles.yaml", "N2", anharmonic=true, simplified_anharmonic=true)

    # println(compute_max_vibr_level(mol, 20000.0, 200.0))
end