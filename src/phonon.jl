#=
  ____   _                                 
 |  _ \ | |__    ___   _ __    ___   _ __  
 | |_) || '_ \  / _ \ | '_ \  / _ \ | '_ \ 
 |  __/ | | | || (_) || | | || (_) || | | |
 |_|    |_| |_| \___/ |_| |_| \___/ |_| |_|
===========================================

This file contains the inner workings for performing calculations with
force constants obtained.
=#

module Phonon

import Logging

include("./parsers.jl")

"""
From this struct we will draw all the data we need to perform a calculation
"""
struct HarmonicScaffolding
    mass_prefactor::Array{Float64}

    # function HarmonicScaffolding(
    #     ebdata::Parsers.ebInputData,
    #     dfptdata::Parsers.dfpt_qeOutputData,
    # ) end
end

"""
    enforce_acoustic_sum_rule!(ifc2_tensor)

Enforce the acoustic sum rule on an tensor of force constants. It can be derived
from Newtons 3rd Law for the atoms in the unitcell at the origin.

The force constants should have the shape `(3, 3, nat, nat, sc[1], sc[2], sc[3])`,
where `nat` are the number of atoms in a unitcell and `sc` is a vector of integers
that correspond to the number of unitcells that make up the supercell over which
the tensor is defined.

# References

  - G. J. Ackland et. al. 1997 "Practical methods in ab initio lattice dynamics"
"""
function enforce_acoustic_sum_rule!(ifc2_tensor::Array{Float64})
    # Grab the number of atoms from the shape
    nat = size(ifc2_tensor)[3]

    for i in 1:3
        for j in 1:3
            for iat in 1:nat
                full_sum = sum(ifc2_tensor[i, j, iat, :, :, :, :])
                ifc2_tensor[i, j, iat, iat, 1, 1, 1] =
                    ifc2_tensor[i, j, iat, iat, 1, 1, 1] - full_sum
            end
        end
    end
    return
end

"""
    build_mass_prefactor(
            numbasisatoms::Int64,
            basisatom2species::Vector{Int64},
            species2mass::Vector{Float64})

Build the mass prefactor that is needed in the computation of the
dynamical matrix.

The mass_prefactor is a matrix of size `(numbasisatoms,  numbasisatoms)`. Its `i`-th diagonal element is just the mass of
`species[basisatom2species[i]]` which we will call `mass[i]` in this
comment. Every other `(i,j)`-th element is `sqrt( mass[i] * mass[j] )`
"""
function build_mass_prefactor(
    numbasisatoms::Int64,
    basisatom2species::Vector{Int64},
    species2mass::Vector{Float64},
)

    # Build a vector that holds the mass for each basisatom
    mass = Vector{Float64}(undef, (numbasisatoms))
    for basisatom in 1:numbasisatoms
        # Grab the species of each basis atom
        basisspecies = basisatom2species[basisatom]
        # Convert the species of each basis atom to its mass
        mass[basisatom] = species2mass[basisspecies]
    end

    mass_prefactor = sqrt.(mass * transpose(mass))

    return mass_prefactor
end

"""
    build_basisconnectors(numbasisatoms::Int64, basis::Matrix{Float64})

Build all vectors that connect each position defined in `basis`.

The tensor that is built is anti-symmetric in the first two indices and
thus contains zero-vectors on its diagonal. The `(i,j)`-th element of
`basisconnectors` contains the difference of `basis[i]` and `basis[j]` with
the `(j,i)`-th element being the negative of the former mentioned.
"""
function build_basisconnectors(numbasisatoms::Int64, basis::Matrix{Float64})

    # Build square (numatoms x numatoms)-matrix of vectors such that every 
    # row consists of all basis vectors of all basisatoms
    basis_dublicate = Array{Float64}(undef, (numbasisatoms, numbasisatoms, 3))

    # Fill the first row
    basis_dublicate[1, :, :] = basis[:, :]
    # Dublicate the first row to the other
    for j in 1:numbasisatoms
        basis_dublicate[j, :, :] = basis_dublicate[1, :, :]
    end

    # Connectors are the difference between the dublicate and the 
    # "transpose" (as for a matrix with vector-elements) of the dublicate 
    basisconnectors = basis_dublicate - permutedims(basis_dublicate, (2, 1, 3))

    return basisconnectors
end

"""
    build_supercell_positions(
        unit_multiplicity_super::Vector{Int64},
        super_multiplicity_ultra::Vector{Int64},
        unit_lattvecs::Matrix{Float64})

Generate an ultracell by repetion of supercells, for which you build
all positions to as well as the positions half length squared.

The supercell is made up of a multiplicity of unitcells in each direction
of `unit_lattvecs`. Analogously the ultracell is made up of a
multiplicity of supercells in each direction of `unit_lattvecs`.

# Important

  - For now `super_multiplicity_ultra` must an even number
  - It is assumed that the `unit_lattvecs` column (2nd index) goes through the
    cartesian coordinates!
  - It is assumed that the orderings in the `unit_multiplicity_super` and
    `unit_lattvecs` align
"""
function build_supercell_points(
    unit_multiplicity_super::Vector{Int64},
    super_multiplicity_ultra::Vector{Int64},
    unit_lattvecs::Matrix{Float64},
)

    # Supercell multiplicity must be even in every direction such that we 
    # have an origin in the middle of the cube
    for i in 1:3
        if isodd(super_multiplicity_ultra[i]) || super_multiplicity_ultra[i] < 2
            Logging.@error "Multiplicity of the supercells in the ultracell is not valid!"
            return
        end
    end

    # Build the lattice vectors of the supercell by stretching the unitcell 
    # lattvecs
    super_lattvecs = Matrix{Float64}(undef, (3, 3))
    for row in 1:3
        super_lattvecs[row, :] = unit_lattvecs[row, :] * unit_multiplicity_super[row]
    end

    # The supercell multiplicity gives the number of cells we will build 
    # inside the ultracell in each lattvec directions. The number of 
    # supercell borders in the i-th direction is 
    # super_multiplicity_ultra[i]+1
    num_supercell_points = prod(super_multiplicity_ultra .+ 1)
    Logging.@info "Phonon.build_supercell_positions: The number of supercell point is" num_supercell_points
    # Saving the supercell positions in the ultracell
    super_points = Matrix{Float64}(undef, (num_supercell_points, 3))
    # Saving the distance of the bisector between every supercell position 
    # and the origin
    super_bisector_dist = Vector{Float64}(undef, num_supercell_points)

    # In Quantum Espresso ultra_range_max is fixed to 2
    ultra_range_max = div.(super_multiplicity_ultra, 2)
    ultra_range = range.(-ultra_range_max, ultra_range_max)
    Logging.@info """
        Phonon.build_supercell_positions: The loop ranges are
    """ super_multiplicity_ultra ultra_range
    # Building the ultracell as a cube with the same number of 
    # supercell_positions in every direction around the origin which is 
    # (0,0,0) here
    j = 1
    for m1 in ultra_range[1]
        for m2 in ultra_range[2]
            for m3 in ultra_range[3]
                # IMPORTANT: in the Quantum Espresso and elphbolt the 
                # origin is skipped at this step here. I want to keep the 
                # origin so it is more verbose (or direct whatever) to skip 
                # it later on and not hidden inside the data structure

                super_points[j, :] = begin
                    super_lattvecs[1, :] * m1 +
                    super_lattvecs[2, :] * m2 +
                    super_lattvecs[3, :] * m3
                end

                super_bisector_dist[j] = begin
                    0.5 * transpose(super_points[j, :]) * super_points[j, :]
                end

                j += 1
            end
        end
    end

    return super_points, super_bisector_dist
end

"""
    build_unitcell_points(
        unit_multiplicity_super::Vector{Int64},
        super_multiplicity_ultra::Vector{Int64},
        unit_lattvecs::Matrix{Float64})

Build a finite lattice of the size of an ultracell specified by the
supercell multiplicity. The lattic is unitcell periodic.

# Important

  - For now `super_multiplicity_ultra` must an even number
  - It is assumed that the `unit_lattvecs` column (2nd index) goes through the
    cartesian coordinates!
  - It is assumed that the orderings in the `unit_multiplicity_super` and
    `unit_lattvecs` align
"""
function build_unitcell_points(
    unit_multiplicity_super::Vector{Int64},
    super_multiplicity_ultra::Vector{Int64},
    unit_lattvecs::Matrix{Float64},
)

    # super_multiplicity_ultra must be confined to even numbers. See the 
    # function Phonon.build_super_points()
    for i in 1:3
        if isodd(super_multiplicity_ultra[i]) || super_multiplicity_ultra[i] < 2
            Logging.@error "Multiplicity of the supercells in the ultracell is not valid!"
            return
        end
    end

    # The number of unitcell points can be calculated as follows. In each 
    # i-th lattvec direction there are super_mult[i]*unit_mult[i]+1 
    # unitcell points because there super_mult[i]*unit_mult[i] cells in 
    # that lattvec direction. For three directions the following gives the 
    # number of unitcell points.
    num_unit_points = begin
        prod(super_multiplicity_ultra .* unit_multiplicity_super .+ 1)
    end
    Logging.@info """
        Phonons.build_unitcell_points: The number of unitcell points is
    """ num_unit_points

    unit_points = Matrix{Float64}(undef, (num_unit_points, 3))

    ultra_range_max = div.(super_multiplicity_ultra, 2) .* unit_multiplicity_super
    ultra_range = range.(-ultra_range_max, ultra_range_max)

    j = 1
    for m1 in ultra_range[1]
        for m2 in ultra_range[2]
            for m3 in ultra_range[3]
                unit_points[j, :] = begin
                    unit_lattvecs[1, :] * m1 +
                    unit_lattvecs[2, :] * m2 +
                    unit_lattvecs[3, :] * m3
                end

                j += 1
            end
        end
    end

    return unit_points
end

"""
    get_origin_ultraconnectors(
        numbasisatoms::Int64,
        unit_points::Matrix{Float64},
        basisconnectors::Array{Float64})

Calculate all vectors from the atoms in the unitcell at (0,0,0) to every
atom in the ultracell (including self-distance).

This is basically building `numbasisatoms` coordinate lists. In each
`i`-th list are the vectors that point to all atoms in the ultracell. But
each `i`-th list is not the same! They differ from each other because
each list is calculated with the origin shifted into the `i`-th atom in the
unitcell at (0,0,0)

# Important

  - It is assumed that unit_points and basisconnectors are given with
    respect to the same (vector space) basis and have the same units
"""
function get_origin_ultraconnectors(
    numbasisatoms::Int64,
    unit_points::Matrix{Float64},
    basisconnectors::Array{Float64},
)

    # Get the number of unit points
    num_unit_points = size(unit_points)[1]

    num_ultraconnectors = numbasisatoms^2 * num_unit_points

    ultraconnectors = Array{Float64}(undef, (num_ultraconnectors, 3))

    s = 1
    for iat in 1:numbasisatoms
        for jat in 1:numbasisatoms
            for k in 1:num_unit_points
                ultraconnectors[s, :] = unit_points[k, :] + basisconnectors[iat, jat, :]

                s += 1
            end
        end
    end

    return ultraconnectors
end

function build_weight_maps()

    # Maybe build the weight maps as if the ultraconnectors really are 
    # `numbasisatoms` lists of coordinates. So one map per list of 
    # coordinates!

end

# End of module Phonon
end
