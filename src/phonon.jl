"""
    <!func signature here!>

# Arguments

  - `lattvecs::Matrix{Float64}`: Collection of direct lattice vectors of the lattice.
    Stored rowwise as in `lattvecs[:,i]` will yield the `i`-th lattice vector.

  - `qpoints_cryst::Matrix{Float64}`: A list of q-points in fractional coordinates for
    which the dynamical matrix will be constructed. Ideally this will be a path based
    on the crystal symmetry. It is assumed that `qpoints_cryst[i,:]` will yield the
    `i`-th q-point in the list.
"""
function build_dynamical_matrix(
    deconvolution::DeconvData,
    qpoints_cryst::Matrix{Float64},
    lattvecs::Matrix{Float64},
    ifc2::Array{Float64,7},
)
    # Extract the data of the deconvolution
    weightmap = deconvolution.weightmap
    unitpoints_qefrac_folded = deconvolution.unitpoints_qefrac_folded
    unitpoints_cart = deconvolution.unitpoints_cart

    # Lattvecs are given in nm we need the reclattvecs in bohr
    # Reciprocal lattice vectors have units of 1/length such that a 
    # multiplication with the Bohr-radius in nm will convert from 1/nm to 1/bohr
    reclattvecs = calc_reciprocal_lattvecs(lattvecs) * a0_nm
    # Change basis of q-points to cartesian in units of Bohr
    qpoints_cart = qpoints_cryst2cart(reclattvecs, qpoints_cryst)

    numqpoints = size(qpoints_cryst, 1)
    numatoms = size(weightmap, 1)
    dynmat = Array{ComplexF64,3}(undef, (numqpoints, numatoms * 3, numatoms * 3))

    # We loop over all atoms in all unitcells in the ultracell. We will check for a 
    # non-vanishing weights at that unitcell
    num_unit_points = size(unitpoints_cart, 1)
    for iat in 1:numatoms
        for jat in 1:numatoms
            for icart in 1:3
                i = mux2to1(iat, icart)
                for jcart in 1:3
                    j = mux2to1(jat, jcart)

                    # dynmat[:,i,j] = calc_fullq_dynmat_element()
                end
            end
        end
    end
end

"""
For all given q-points calculate one term of the lattice sum for one submatrix of the
dynamical matrix.

For one q-point the dynamical tensor is a rank-4 tensor indexed by two atom indeces
τ,τ' (`iat` and `jat`) and two cartesian indeces α, α' (`ipol` and `jpol`). Solving
for frequencies and vibration eigenvectors in one diagonalization per q-point
requires to mux the indices τ and α into one index `idim` and τ' and α' analogously
into `jdim` using `mux2to1()`. For a given q-point this converts the dynamical tensor
(rank-4) into the dynamical matrix (rank-2). If we pick a q-point and fix τ and τ' we
get a submatrix of the dynamical matrix.

# Arguments

  - `qpoints_cart::Matrix{Float64}`: A list of q-points in cartesian coordinates.
    Ideally this will be a path based on the crystal symmetry. It is assumed that
    `qpoints_cart[i,:]` will yield the `i`-th q-point in the list.
"""
function calc_fullq_dynmat_element(i::Int64, j::Int64, qpoints_cart::Matrix{Float64})
end

"""
    mux2to1(i::Int64, j::Int64)

Mux two indeces into one.

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
    qpoints_cryst::Matrix{Float64},
    reclattvecs::Matrix{Float64},
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
