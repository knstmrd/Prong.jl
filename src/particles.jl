struct Atom
    mass::Float64
    formation_energy::Float64
    degeneracy::UInt8
end

struct Molecule
    mass::Float64
    dissociation_energy::Float64
    formation_energy::Float64
    vibr_energy0::Float64

    rotational_symmetry::UInt8
    
    # infinite_harmonic::Bool
    anharmonic::Bool
    homonuclear::Bool
    use_Treanor::Bool
    continue_Treanor_with_Boltzmann::Bool  # if Tv>T, do we continue the Treanor distribution with a Boltzmann one?

    degeneracy::UInt8

    n_vibr::UInt16
    n_rot::UInt16

    vibrational_levels::Array{Float64,1}
    vibrational_levels_mult1::Array{Float64,1}
    vibrational_energy::Array{Float64,1}

    rotational_energy::Array{Float64,1}
    rotational_degeneracies::Array{Float64,1}
end

function create_atom(filename::String; include_electronic_degeneracy::Bool=true,
                     formation_energy::String="")

end

function create_molecule(filename::String; anharmonic::Bool=true,
                         use_Treanor::Bool=true, continue_Treanor_with_Boltzmann::Bool=true,
                         include_electronic_degeneracy::Bool=true,
                         formation_energy::String="")
    if include_electronic_degeneracy == false
        degeneracy = 1
    end
end