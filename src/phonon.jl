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
function enforce_acoustic_sum_rule!(ifc2_tensor::Array{Float64,7})
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
