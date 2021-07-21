using YAML

struct Atom
    name::String
    mass::Float64
    formation_energy::Float64
    degeneracy::UInt8
    gupta_yos_coefficients::Array{Float64,2}
end

struct Molecule
    name::String
    mass::Float64
    dissociation_energy::Float64
    formation_energy::Float64
    vibr_energy0::Float64
    vibr_energy1::Float64

    rotational_symmetry::UInt8
    
    # infinite_harmonic::Bool
    anharmonic::Bool
    use_Treanor::Bool
    continue_Treanor_with_Boltzmann::Bool  # if Tv>T, do we continue the Treanor distribution with a Boltzmann one?
    frequency::Float64
    anharmonic_ratio::Float64

    degeneracy::UInt8

    n_vibr::UInt16  # number of elements in vibrational energy array = number of max vibrational level + 1
    n_rot::UInt16

    vibrational_levels::Array{Float64,1}
    vibrational_levels_mult1::Array{Float64,1}
    vibrational_levels_mult_energy::Array{Float64,1}
    vibrational_energy::Array{Float64,1}

    rotational_energy::Array{Float64,1}
    rotational_degeneracies::Array{Float64,1}
    gupta_yos_coefficients::Array{Float64,2}
end

function create_atom(filename::String, name::String; include_electronic_degeneracy::Bool=true)

    data = YAML.load(open(filename))
    mass = data[name]["Mass, kg"]

    degeneracy = data[name]["Statistical weight"][1]
    if include_electronic_degeneracy == false
        degeneracy = 1
    end

    fe = data[name]["Formation energy, J"]

    gupta_yos_coefficients = reshape(data[name]["Gupta Yos thermodynamic curve fit coefficients"],(7,5))

    return Atom(name, mass, fe, degeneracy, gupta_yos_coefficients)
end

function create_molecule(filename::String, name::String; anharmonic::Bool=true, simplified_anharmonic::Bool=true,
                         use_Treanor::Bool=true, continue_Treanor_with_Boltzmann::Bool=true,
                         include_electronic_degeneracy::Bool=true)
    data = YAML.load(open(filename))
    mass = data[name]["Mass, kg"]

    degeneracy = data[name]["Statistical weight"][1]

    if include_electronic_degeneracy == false
        degeneracy = 1
    end

    fe = data[name]["Formation energy, J"]

    dissociation_energy = data[name]["Dissociation energy, J"][1]

    rotational_symmetry = data[name]["Factor of symmetry"]

    we = data[name]["Frequency of vibrations (we), m^-1"][1]

    wexe = 0.0
    weye = 0.0
    weze = 0.0
    anharmonic_ratio = 0.0

    if anharmonic
        wexe = data[name]["wexe, m^-1"][1]
        anharmonic_ratio = wexe / we

        if !simplified_anharmonic
            weye = data[name]["weye, m^-1"][1]
            weze = data[name]["weye, m^-1"][1]
        end
    end

    i = 0
    tmp = 0.0

    ve_arr = []
    vl_arr = []

    while (tmp < dissociation_energy)
        tmp = constants.h * constants.c * (we * (i + 0.5) - wexe * (i + 0.5)^2 
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

    rot_be = data[name]["Be, m^-1"][1]

    j = 0
    tmp = 0.0

    re_arr = []
    rd_arr = []

    while (tmp < dissociation_energy)
        tmp = constants.h * constants.c * (rot_be * j * (j+1));
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

    use_Tr = use_Treanor
    cont_Tr_B = continue_Treanor_with_Boltzmann
    if !anharmonic
        use_Tr = false
    end

    if !use_Tr
        cont_Tr_B = false
    end

    gupta_yos_coefficients = reshape(data[name]["Gupta Yos thermodynamic curve fit coefficients"],(7,5))
    
    return Molecule(name, mass, dissociation_energy, fe, ve0, ve_arr[2], rotational_symmetry, anharmonic, use_Tr, cont_Tr_B,
                    we, anharmonic_ratio, degeneracy, n_vibr, n_rot, vl_arr, vl_x1_arr, vl_xve_arr, ve_arr, re_arr, rd_arr,
                    gupta_yos_coefficients)
end