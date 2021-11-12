using YAML
using Unitful
using PhysicalConstants.CODATA2018

struct Atom
    name::String
    mass::typeof(1.0u"kg")
    formation_energy::typeof(1.0u"J")
    degeneracy::UInt8
    gupta_yos_coefficients::Array{Float64,2}
end

struct Molecule
    name::String
    mass::typeof(1.0u"kg")
    dissociation_energy::typeof(1.0u"J")
    formation_energy::typeof(1.0u"J")
    vibr_energy0::typeof(1.0u"J")
    vibr_energy1::typeof(1.0u"J")

    rotational_symmetry::UInt8
    
    # infinite_harmonic::Bool
    anharmonic::Bool
    # use_Treanor::Bool
    # continue_Treanor_with_Boltzmann::Bool  # if Tv>T, do we continue the Treanor distribution with a Boltzmann one?
    frequency::typeof(1.0u"m^-1")
    anharmonic_ratio::Float64

    degeneracy::UInt8

    n_vibr::UInt16  # number of elements in vibrational energy array = number of max vibrational level + 1
    n_rot::UInt16

    vibrational_levels::Array{Float64,1}
    vibrational_levels_mult1::Array{typeof(1.0u"J"),1}  # i * vibr_energy(first level)
    vibrational_levels_mult_energy::Array{typeof(1.0u"J"),1}  # i * vibr_energy(level i)
    Δvibrational_anharmonic::Array{typeof(1.0u"J"),1}
    vibrational_energy::Array{typeof(1.0u"J"),1}

    rotational_energy::Array{typeof(1.0u"J"),1}
    rotational_degeneracies::Array{Float64,1}
    gupta_yos_coefficients::Array{Float64,2}
end

struct VibrationalDistribution
    use_Treanor::Bool
    continue_Treanor_with_Boltzmann::Bool
end

function create_atom(filename::String, name::String; include_electronic_degeneracy::Bool=true)

    data = YAML.load(open(filename))
    mass = data[name]["Mass, kg"] * 1.0u"kg"

    degeneracy = data[name]["Statistical weight"][1]
    if include_electronic_degeneracy == false
        degeneracy = 1
    end

    fe = data[name]["Formation energy, J"] * 1.0u"J"

    gupta_yos_coefficients = reshape(data[name]["Gupta Yos thermodynamic curve fit coefficients"],(7,5))

    return Atom(name, mass, fe, degeneracy, gupta_yos_coefficients)
end

function create_molecule(filename::String, name::String; anharmonic::Bool=true, simplified_anharmonic::Bool=true,
                         include_electronic_degeneracy::Bool=true)
    data = YAML.load(open(filename))
    mass = data[name]["Mass, kg"] * 1.0u"kg"

    degeneracy = data[name]["Statistical weight"][1]

    if include_electronic_degeneracy == false
        degeneracy = 1
    end

    fe = data[name]["Formation energy, J"] * 1.0u"J"

    dissociation_energy = data[name]["Dissociation energy, J"][1] * 1.0u"J"

    rotational_symmetry = data[name]["Factor of symmetry"]

    we = data[name]["Frequency of vibrations (we), m^-1"][1] * 1.0u"m^-1"

    wexe = 0.0u"m^-1"
    weye = 0.0u"m^-1"
    weze = 0.0u"m^-1"
    anharmonic_ratio = 0.0

    if anharmonic
        wexe = data[name]["wexe, m^-1"][1] * 1.0u"m^-1"
        anharmonic_ratio = wexe / we

        if !simplified_anharmonic
            weye = data[name]["weye, m^-1"][1] * 1.0u"m^-1"
            weze = data[name]["weye, m^-1"][1] * 1.0u"m^-1"
        end
    end

    i = 0
    tmp = 0.0u"J"

    ve_arr = []
    vl_arr = []

    while (tmp < dissociation_energy)
        tmp = PlanckConstant * SpeedOfLightInVacuum	* (we * (i + 0.5) - wexe * (i + 0.5)^2 
                                                          + weye * (i + 0.5)^3
                                                          + weze * (i + 0.5)^4)
        if (tmp < dissociation_energy)
            if i == 0
                push!(ve_arr, tmp)
                push!(vl_arr, i)
                i += 1
            else
                if tmp > ve_arr[i]
                    push!(ve_arr, tmp)
                    push!(vl_arr, i)
                    i += 1
                else
                    println("Non-monotonous vibrational energy spectrum computed for ", name, ", ", tmp, ", ", ve_arr[i])
                    break
                end
            end
        end
    end

    n_vibr = length(ve_arr)
    ve0 = ve_arr[1]
    ve_arr .-= ve0

    vl_x1_arr = vl_arr .* ve_arr[2]
    vl_xve_arr = vl_arr .* ve_arr

    rot_be = data[name]["Be, m^-1"][1] * 1.0u"m^-1"

    j = 0
    tmp = 0.0u"J"

    re_arr = []
    rd_arr = []

    while (tmp < dissociation_energy)
        tmp = PlanckConstant * SpeedOfLightInVacuum * (rot_be * j * (j+1));
        if (tmp < dissociation_energy)
            if j == 0
                push!(re_arr, tmp)
                push!(rd_arr, 2 * j + 1)
                j += 1
            else
                if tmp > re_arr[j]
                    push!(re_arr, tmp)
                    push!(rd_arr, 2 * j + 1)
                    j += 1
                else
                    break
                end
            end
        end
    end

    n_rot = length(re_arr)

    gupta_yos_coefficients = reshape(data[name]["Gupta Yos thermodynamic curve fit coefficients"],(7,5))

    Δvibrational_anharmonic = ve_arr - vl_x1_arr
    # mol.vibrational_energy .- mol.vibrational_levels_mult1
    
    return Molecule(name, mass, dissociation_energy, fe, ve0, ve_arr[2], rotational_symmetry, anharmonic,
                    we, anharmonic_ratio, degeneracy, n_vibr, n_rot, vl_arr, vl_x1_arr, vl_xve_arr, Δvibrational_anharmonic,
                    ve_arr, re_arr, rd_arr,
                    gupta_yos_coefficients)
end

function create_vibrational_distribution(anharmonic::Bool, use_Treanor::Bool, continue_Treanor_with_Boltzmann::Bool)
    use_Tr = use_Treanor
    cont_Tr_B = continue_Treanor_with_Boltzmann
    
    if !anharmonic
        use_Tr = false
    end

    if !use_Tr
        cont_Tr_B = false
    end

    return VibrationalDistribution(use_Tr, cont_Tr_B)
end