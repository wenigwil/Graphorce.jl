struct LatticeVibrations
    fullq_freqs::Matrix{Float64}
    qpoints_cryst::Matrix{Float64}
    stitches::BitVector

    function LatticeVibrations(
        ebinputfilepath::AbstractString,
        qeifc2outputfilepath::AbstractString,
        highsympoints_cryst...;
        numpoints_per_section::Int64 = 50,
    )
        # Read necessary data from file for computation
        ebdata = ebInputData(ebinputfilepath)
        qedata = qeIfc2Output(qeifc2outputfilepath)
        qpoints_cryst, stitches = points_to_path(
            highsympoints_cryst...;
            numpoints_per_section = numpoints_per_section,
            return_stitches = true,
        )
        # Create a deconvolution instance
        deconvolution = DeconvData(ebdata)

        lattvecs = ebdata.crystal_info["lattvecs"]
        basisatoms2species = ebdata.crystal_info["atomtypes"]
        species2masses = qedata.properties["species2masses"]
        ifc2 = qedata.properties["ifc2"]
        weightmap = deconvolution.weightmap
        uqf = deconvolution.unitpoints_qefrac_folded
        unitpoints_cart = deconvolution.unitpoints_cart

        numqpoints = size(qpoints_cryst, 1)
        numatoms = size(basisatoms2species, 1)

        enforce_acoustic_sum_rule!(ifc2)

        fullq_freqs = Matrix{Float64}(undef, (numqpoints, 3 * numatoms))
        for iq in 1:numqpoints
            dynmat = build_dynamical_matrix(
                weightmap,
                uqf,
                unitpoints_cart,
                lattvecs,
                ifc2,
                basisatoms2species,
                species2masses,
                qpoints_cryst[iq, :],
            )

            force_hermiticity!(dynmat)

            dynmat = LinAlg.Hermitian(dynmat)

            eigvals = LinAlg.eigvals(dynmat)

            fullq_freqs[iq, :] = copysign.(sqrt.(abs.(eigvals)), eigvals)

            if iszero(qpoints_cryst[iq, :])
                fullq_freqs[iq, 1:3] .= [0.0, 0.0, 0.0]
            end
        end

        # Unit conversion
        fullq_freqs .*= RydtoTHz

        new(fullq_freqs, qpoints_cryst, stitches)
    end
end

function force_hermiticity!(mat::Matrix{ComplexF64})
    mat .= 1 / 2 * (mat + transpose(conj(mat)))
end

"""
The dynamical tensor is a rank-5 tensor with two atomic indices τ, τ' (inclusive
range 1 to `numatoms`), two cartesian indices α, α' (inclusive range 1 to 3) and one
q-point index `iq` (inclusive range 1 to `numqpoints`). The tensor will be
represented as a rank-3 tensor by muxing the two index-pairs τ, α and τ', α' into two
indices `i` and `j`. For a specified `iq` we then get a square
`3*numatoms`x`3*numatoms`-matrix.
"""
function build_dynamical_matrix(
    weightmap::Array{Float64,3},
    uqf::Array{Int64,4},
    unitpoints_cart::Matrix{Float64},
    lattvecs::Matrix{Float64},
    ifc2::Array{Float64,7},
    basisatom2species::Vector{Int64},
    species2mass::Vector{Float64},
    qpoint_cryst::Vector{Float64},
)
    numatoms = size(basisatom2species, 1)
    dynmat = zeros(ComplexF64, (3 * numatoms, 3 * numatoms))
    numunitpoints = size(unitpoints_cart, 1)

    # Lattvecs are given in nm we need the reclattvecs in bohr
    # Reciprocal lattice vectors have units of 1/length such that a 
    # multiplication with the Bohr-radius in nm will convert from 1/nm to 1/bohr
    reclattvecs = calc_reciprocal_lattvecs(lattvecs) * a0_nm

    qpoint_cart = reclattvecs * qpoint_cryst

    mass_prefactor = build_mass_prefactor(basisatom2species, species2mass)
    for iat in 1:numatoms
        for jat in 1:numatoms
            for icart in 1:3
                i = mux2to1(iat, icart)
                for jcart in 1:3
                    j = mux2to1(jat, jcart)

                    for l in 1:numunitpoints
                        if weightmap[l, jat, iat] > 0
                            dynmat[i, j] += begin
                                ifc2[
                                    icart,
                                    jcart,
                                    iat,
                                    jat,
                                    uqf[l, jat, iat, 1],
                                    uqf[l, jat, iat, 2],
                                    uqf[l, jat, iat, 3],
                                ] *
                                exp(
                                    im * LinAlg.dot(
                                        qpoint_cart,
                                        unitpoints_cart[l, :],
                                    ),
                                ) *
                                weightmap[l, jat, iat]
                            end
                        end
                    end

                    dynmat[i, j] /= mass_prefactor[iat, jat]
                end
            end
        end
    end

    return dynmat
end

"""
    mux2to1(i::Int64, j::Int64)

Mux two indeces into one. This is only a valid muxing method for 1-based indices.

A very old man told me, it was a good idea.
"""
function mux2to1(i::Int64, j::Int64)
    return (i - 1) * 3 + j
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
  - `reclattvecs::Matrix{Float64}`: Reciprocal lattice vectors. `reclattvecs[:,i]`
    will yield the `i`-th reciprocal lattice vector.
"""
function qpoints_cryst2cart(
    reclattvecs::Matrix{Float64},
    qpoints_cryst::Matrix{Float64},
)
    numqpoints = size(qpoints_cryst, 1)
    qpoints_cart = zeros(Float64, (numqpoints, 3))
    for iq in 1:numqpoints
        qpoints_cart[iq, :] = reclattvecs * qpoints_cryst[iq, :]
    end

    return qpoints_cart
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
    nat = size(ifc2_tensor, 3)

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

The mass_prefactor is a matrix of size `(numbasisatoms,  numbasisatoms)`. Its `i`-th
diagonal element is just the mass of `species[basisatom2species[i]]` which we will
call `mass[i]` in this comment. Every other `(i,j)`-th element is
`sqrt( mass[i] * mass[j] )`
"""
function build_mass_prefactor(
    basisatom2species::Vector{Int64},
    species2mass::Vector{Float64},
)
    numbasisatoms = size(basisatom2species, 1)
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
