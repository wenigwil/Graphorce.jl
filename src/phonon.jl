struct DynamicalMatrix
    qpoints_cart::Vector{Float64}
    mat::Array{ComplexF64,3}

    function DynamicalMatrix(
        qpoints_cryst::Vector{Float64},
        lattvecs::Matrix{Float64},
    )

        # Lattvecs are given in nm we need the reclattvecs in bohr
        # Reciprocal lattice vectors have units of 1/length such that a multiplication 
        # with the Bohr-radius in nm will convert from 1/nm to 1/bohr
        reclattvecs = transpose(calc_reciprocal_lattvecs(lattvecs) * a0_nm)




        new()
    end
end

"""
    qpoints_cryst2cart(
        qpoints_cryst::Matrix{Float64},
        reclattvecs::Matrix{Float64})

Convert a list of reciprocal lattice vectors given in crystal coordinates to ones
that are given in cartesian coordinates with units of `reclattvecs`.

# Arguments

  - `qpoints_cryst::Matrix{Float64}`: A list of q-points given in crystal
    coordinates. `qpoints_cryst[i,:]` should yield the `i`-th q-point in the list.
  - `reclattvecs::Matrix{Float64}`: Reciprocal lattice vectors. `reclattvecs[i,:]`
    will yield the `i`-th reciprocal lattice vector.
"""
function qpoints_cryst2cart(
    qpoints_cryst::Matrix{Float64},
    reclattvecs::Matrix{Float64},
)
    numqpoints = size(qpoints_cryst)[1]
    qpoints_cart = zeros(Float64, (numqpoints, 3))
    for iq in 1:numqpoints
        qpoints_cart[iq, :] = reclattvecs * qpoints_cryst[iq, :]
    end
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
