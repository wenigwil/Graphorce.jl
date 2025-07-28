struct DeconvUtils
    weightmap::Array{Float64,3}
    ifc2_2nd_atom_unitaddr::Array{Int64,4}

    function DeconvUtils(ebdata::ebInputData; super_multiplicity_ultra = [4, 4, 4])
        numbasisatoms = ebdata.allocations["numatoms"]
        # numbasisatomsx3 basis in fractional coordiantes
        basis_frac = transpose(ebdata.crystal_info["basis"])
        # lattice vectors given as transpose([a1 a2 a3]) in Bohr
        lattvecs = transpose(ebdata.crystal_info["lattvecs"]) / a0_nm
        # Get the basis in cartesian coords and convert to Bohr
        # Will give a numbasisatomsx3 matrix
        basis_cart = (basis_frac * lattvecs)
        basiscons_cart = build_basisconnectors(numbasisatoms, basis_cart)

        unit_multiplicity_super = ebdata.numerics["qmesh"]
        # Get the unitcell points in cartesian coords in units of Bohr
        unit_points_cart = build_unitcell_points(
            unit_multiplicity_super,
            super_multiplicity_ultra,
            lattvecs,
        )

        # unit_points and basiscons must be in the same space and units
        shiftercons =
            get_shiftercons(numbasisatoms, unit_points_cart, basiscons_cart)
        # super_sqmods stands for the squared modulus of super_points
        super_points, super_sqmods = build_supercell_points(
            unit_multiplicity_super,
            super_multiplicity_ultra,
            lattvecs,
        )

        weightmap = get_weight_map(shiftercons, super_points, super_sqmods)
        unit_points_frac = build_unitcell_points_frac(
            unit_multiplicity_super,
            super_multiplicity_ultra,
        )

        ifc2_2ndatom_unitaddr = get_ifc2_2nd_atom_unitcell_addr(
            weightmap,
            unit_points_frac,
            unit_multiplicity_super,
        )

        new(weightmap, ifc2_2ndatom_unitaddr)
    end
end

"""
    build_basisconnectors(numbasisatoms::Int64, basis::Matrix{Float64})

Build all vectors that connect each position defined in `basis`.

The tensor that is built is anti-symmetric in the first two indices and
thus contains zero-vectors on its diagonal. The `(i,j)`-th element of
`basisconnectors` contains the difference of `basis[i]` and `basis[j]` with
the `(j,i)`-th element being the negative of the former mentioned.

# Arguments

  - `numbasisatoms::Int64`: Number of atoms in the crystal basis
  - `basis::Matrix{Float64}`: Collection of basis vectors of atoms in the
    unitcell. Size should be (`numbasis`, 3) with the basis vectors forming
    the rows of `basis`.

# Output

  - `basisconnectors::Array{Float64,3}`: Element `[i,j,:]` will yield the vector
    from atom j to atom i, where both atoms are inside of the same unitcell.
"""
function build_basisconnectors(numbasisatoms::Int64, basis::Matrix{Float64})

    # Build square (numatoms x numatoms)-matrix of vectors such that every 
    # row consists of all basis vectors of all basisatoms
    basis_dublicate = Array{Float64,3}(undef, (numbasisatoms, numbasisatoms, 3))

    basis_dublicate[1, :, :] = basis[:, :]
    # Dublicate the first row to the other rows
    for j in 1:numbasisatoms
        basis_dublicate[j, :, :] = basis_dublicate[1, :, :]
    end

    # Connectors are the difference between the dublicate and the 
    # "transpose" (as for a matrix with vector-elements) of the dublicate
    basisconnectors = permutedims(basis_dublicate, (2, 1, 3)) - basis_dublicate

    return basisconnectors
end

"""
    build_supercell_points(
        unit_multiplicity_super::Vector{Int64},
        super_multiplicity_ultra::Vector{Int64},
        unit_lattvecs::Matrix{Float64})

Generate an ultracell by repetion of supercells, for which you build
all positions to as well as the positions half length squared.

The supercell is made up of a multiplicity of unitcells in each direction
of `unit_lattvecs` with the origin at the a corner of the supercell. The
ultracell has its origin at the center.

# Arguments

  - `unit_multiplicity_super::Vector{Int64}`: Multiplicities of unitcells inside
    of a supercell in each of the lattice vector directions. The `j`-th element
    corresponds to the multiplicity in the lattice vector direction given by
    `unit_lattvecs[j,:]`.
  - `super_multiplicity_ultra::Vector{Int64}`: `super_multiplicity_ultra[j]`
    corresponds to the number of supercells in the direction of
    `unit_lattvecs[j,:]` and `-unit_lattvecs[j,:]` so the number of supercells
    along that direction is `2*super_multiplicity_ultra[j]`. Every element shall
    be even and not smaller than 2.
  - `unit_lattvecs::Matrix{Float64}`: Three dimensional lattice vectors. The
    `j`-th lattice vector corresponds to `unit_lattvecs[j,:]`.

# Output

  - `super_points::Matrix{Float64}`: List of the origins of all supercells
    inside the ultracell
  - `super_sqmods::Vector{Float64}`: Squared modulus of every vector in
    `super_points`
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
    Logging.@debug """build_supercell_positions: The number of supercell point is
    """ num_supercell_points
    # Saving the supercell positions in the ultracell
    super_points = Matrix{Float64}(undef, (num_supercell_points, 3))
    # Saving the SQUARED distance between supercell positions and origin
    super_point_sqmods = Vector{Float64}(undef, num_supercell_points)

    # In Quantum Espresso ultra_range_max is fixed to 2
    ultra_range_max = div.(super_multiplicity_ultra, 2)
    ultra_range = range.(-ultra_range_max, ultra_range_max)
    Logging.@debug """
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
                # origin so it is more verbose (or direct, what have you) 
                # to skip it later on and not hidden inside the data 
                # structure

                super_points[j, :] = begin
                    super_lattvecs[1, :] * m1 +
                    super_lattvecs[2, :] * m2 +
                    super_lattvecs[3, :] * m3
                end

                super_point_sqmods[j] =
                    transpose(super_points[j, :]) * super_points[j, :]

                j += 1
            end
        end
    end

    return super_points, super_point_sqmods
end

"""
    build_unitcell_points(
        unit_multiplicity_super::Vector{Int64},
        super_multiplicity_ultra::Vector{Int64},
        unit_lattvecs::Matrix{Float64})

Build a finite lattice of the size of an ultracell specified by the
supercell multiplicity. The lattice is unitcell periodic.

# Arguments

  - `unit_multiplicity_super::Vector{Int64}`: Multiplicities of unitcells inside
    of a supercell in each of the lattice vector directions. The `j`-th element
    corresponds to the multiplicity in the lattice vector direction given by
    `unit_lattvecs[j,:]`.
  - `super_multiplicity_ultra::Vector{Int64}`: Multiplicities of supercells inside
    of the ultracell in positive lattice vector directions. Should coincide with
    the multiplicities supplied in `build_supercell_points()`. Every element must
    be even and may not be smaller than 2.
  - `unit_lattvecs::Matrix{Float64}`: Three dimensional lattice vectors. The
    `j`-th lattice vector corresponds to `unit_lattvecs[j,:]`.

# Output

  - `unit_points::Matrix{Float64}`: A list of all unitcell points inside the
    ultracell that is build up of supercells from `build_supercell_points()`. The
    points are given in the same units as `unit_lattvecs`.
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
    Logging.@debug """
        build_unitcell_points: The number of unitcell points is
    """ num_unit_points

    unit_points = Matrix{Float64}(undef, (num_unit_points, 3))

    # Same concept from Phonon.build_supercell_points()
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
    get_shiftercons(
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

# Arguments

  - `numbasisatoms::Int64`: Number of atoms in one unitcell
  - `unit_points::Matrix{Float64}`: Output from `build_unitcell_points()`.
    Should be in same vector space basis and units as `basisconnectors`
  - `basisconnectors::Array{Float,3}`: Output from `build_basisconnectors()`.
    Should be in the same vector space basis and units as `unit_points`

# Output

  - `shiftercons::Array{Float64,4}`: `shiftercons[k,j,i,:]` will yield the
    vector pointing from atom `i` out of the origin unitcell to atom `j` in
    unitcell `k`. `k` runs through all unitcells in the ultracell but it is a
    muxed index that can be demuxed with the loops from `build_unitcell_points()`.
"""
function get_shiftercons(
    numbasisatoms::Int64,
    unit_points::Matrix{Float64},
    basisconnectors::Array{Float64,3},
)

    # Get the number of unit points
    num_unit_points = size(unit_points)[1]

    # The number of shiftercons is simple to calculate if you think about 
    # how many elements there are in the above mentioned coordinate 
    # lists. We shift one atom from the (0,0,0)-unitcell into the 
    # origin and calculate all vectors pointing from this new origin to 
    # every other atom in the ultracell. So each coordinate list has
    # `numbasisatoms * num_unit_points` elements and there are 
    # `numbasisatoms` coordinate list.
    # num_shiftercons = numbasisatoms^2 * num_unit_points
    shiftercons =
        Array{Float64,4}(undef, (num_unit_points, numbasisatoms, numbasisatoms, 3))

    # Put atom iat from (0,0,0) into the origin
    for iat in 1:numbasisatoms
        for jat in 1:numbasisatoms
            for unit_addr in 1:num_unit_points
                # Compute the vector from atom iat to atom jat in some 
                # unitcell at unit_points[unit_addr]
                shiftercons[unit_addr, jat, iat, :] =
                    unit_points[unit_addr, :] + basisconnectors[jat, iat, :]
            end
        end
    end

    return shiftercons
end

"""
    calc_parlinski_weight(
        shiftercon::Vector{Float64},
        super_points::Matrix{Float64},
        super_point_sqmods::Vector{Float64};
        epsilon::Float64 = 1.0e-6)

Calculate the weight associated to a given shiftercon.

`shiftercon` points from the origin of the ultracell to some point P.
All bisector planes that can be generated from the `super_points` are
checked against a condition for P. That is whether point P is located
in the half space that contains the origin (inside), in the other
half space (outside) or if it is an element of the dividing plane.

  - `weight=0` if point P is considered outside ANY bisector plane
  - `weight=1` if point P is considered inside ALL bisector planes
  - `weight=1/(degen+1)` if P is considered as an element of `degen` planes

In the third case `1/weight` corresponds to the number of
translationally equivalent points to P. The restrictions ANY and ALL in
cases one and two cause to restrict the whole treatment to a Wigner-Seitz
cell around the origin built with the `super_points`.

# Arguments

  - `shiftercon::Vector{Float64}`: One element of `get_shiftercons()` which
    is specified by an index like `shiftercons[k,j,i,:]`.
  - `super_points::Matrix{Float64}`: See output from `build_supercell_points()`
  - `super_sqmods::Vector{Float64}`: See output from `build_supercell_points()`
  - `epsilon::Float64 = 1.0e-6`: Small value for checking whether a point is
    element of a dividing plane or not

# References

  - K. Parlinski 1999 "Calculation of the Phonon Dispersion Curves by the Direct Method"
  - K. Parlinski et. al. 1997 "First-Principles Determination of the Soft Mode in Cubic ZrO2"
"""
function calc_parlinski_weight(
    shiftercon::Vector{Float64},
    super_points::Matrix{Float64},
    super_point_sqmods::Vector{Float64};
    epsilon::Float64 = 1.0e-6,
)::Float64

    # Get the number of super_points there are 
    num_super_points = size(super_points)[1]

    weight = 0.0
    degen = 1
    for super_addr in 1:num_super_points
        # The origin will always fulfill the third case from above so we 
        # skip it here. This is also mentioned in build_supercell_points()
        if iszero(super_points[super_addr, :])
            continue
        end

        # Define numerical value that is checked against for whether the 
        # shiftercon is element of a plane or else.
        check_on_plane = begin
            transpose(shiftercon) * super_points[super_addr, :] -
            0.5 * super_point_sqmods[super_addr]
        end

        # Check whether shiftercon points outside of ANY bisector plane. If 
        # so we immediatly return with weight=0
        # This discards all weight calculations beyond the origin-
        # Wigner-Seitz cell of the superpoints  
        if check_on_plane > epsilon
            return weight
        end

        # Check whether shiftercon is element of a plane.
        if abs(check_on_plane) < epsilon
            # println(
            #     shiftercon,
            #     " is element of plane-generating point ",
            #     super_points[super_addr, :],
            # )
            degen += 1
        end

        # Checking against whether shiftercon is considered inside of ALL 
        # bisector half-spaces is not necessary as this is contained in the 
        # case for no incrementation in degen
    end

    # This stills return 1 if given shiftercon is inside ALL half spaces
    weight = 1 / degen
    return weight
end

"""
    get_weight_map(
        shiftercons::Array{Float64,4},
        super_points::Matrix{Float64},
        super_point_sqmods::Vector{Float64},
    )

Calculate weights with the `calc_parlinski_weight`-function for all `shiftercons` built by
the `get_shiftercons`-function and save it into an array.
"""
function get_weight_map(
    shiftercons::Array{Float64,4},
    super_points::Matrix{Float64},
    super_point_sqmods::Vector{Float64},
)

    # Get the number of unit_points
    num_unit_points, numbasisatoms = size(shiftercons)[1:2]

    weight_map =
        Array{Float64,3}(undef, (num_unit_points, numbasisatoms, numbasisatoms))

    # It is important that this order of loops matches the ordering of the 
    # ones for the shiftercons' creation.
    for iat in 1:numbasisatoms
        for jat in 1:numbasisatoms
            for unit_addr in 1:num_unit_points
                weight_map[unit_addr, jat, iat] = calc_parlinski_weight(
                    shiftercons[unit_addr, jat, iat, :],
                    super_points,
                    super_point_sqmods,
                )
            end
        end
    end

    return weight_map
end

"""
    build_unitcell_points_frac(
        super_multiplicity_ultra::Vector{Int64},
        unit_multiplicity_super::Vector{Int64})

Equivalent to `build_unitcell_points()` yet in fractional coordinates.
Arguments and output are analogous.
"""
function build_unitcell_points_frac(
    super_multiplicity_ultra::Vector{Int64},
    unit_multiplicity_super::Vector{Int64},
)::Matrix{Int64}
    num_unit_points = begin
        prod(super_multiplicity_ultra .* unit_multiplicity_super .+ 1)
    end

    unit_points_frac = Matrix{Int64}(undef, (num_unit_points, 3))
    ultra_range_max = div.(super_multiplicity_ultra, 2) .* unit_multiplicity_super
    ultra_range = range.(-ultra_range_max, ultra_range_max)

    j = 1
    for m1 in ultra_range[1]
        for m2 in ultra_range[2]
            for m3 in ultra_range[3]
                unit_points_frac[j, :] = [m1, m2, m3]
                j += 1
            end
        end
    end

    return unit_points_frac
end

"""
    get_ifc2_2nd_atom_unitcell_addr(
        weight_map::Array{Float64,3},
        unit_points_frac::Matrix{Int64},
        unit_multiplicity_super::Vector{Int64})

Fold all vectors in `unit_points_frac` back into the supercell in the origin and
convert them to a 1-based numbering using `fold_to_qe_origin()`.

With every shiftercon (see explanation a `get_shiftercons()`), for which a weight>=0
was computed, we need to capture the unitcell point in fractional coordinates of the
unshifted atom and fold that back into the origin-supercell. The folding is only
happening if the non-shifted atom is located within a unitcell that is beyond the
origin-supercell.
"""
function get_ifc2_2nd_atom_unitcell_addr(
    weight_map::Array{Float64,3},
    unit_points_frac::Matrix{Int64},
    unit_multiplicity_super::Vector{Int64},
)
    num_unit_points, numbasisatoms = size(weight_map)[1:2]
    unitfold_map = zeros(Int64, (num_unit_points, numbasisatoms, numbasisatoms, 3))

    for iat in 1:numbasisatoms
        for jat in 1:numbasisatoms
            for unit_addr in 1:num_unit_points
                if weight_map[unit_addr, jat, iat] > 0.0
                    unitfold_map[unit_addr, jat, iat, :] = fold_to_qe_origin(
                        unit_points_frac[unit_addr, :],
                        unit_multiplicity_super,
                    )
                end
            end
        end
    end

    return unitfold_map
end

"""
    fold_to_qe_origin(
        unit_point_frac::Vector{Int64},
        unit_multiplicity_super::Vector{Int64})

Fold a unitcell point in fractional coordinates back into a unitcell point that
is inside the supercell at the origin whilst converting to a 1-based (QE-like)
coordinate system inside this supercell.

If the unitcell-multiplicity in one direction was 6, then this
function will map the unit point in that direction NOT onto the
integer interval `[0,5]` but `[1,6]`!

# Arguments

  - `unit_point_frac::Vector{Int64}`: Output from `build_unitcell_points_frac()`
  - `unit_multiplicity_super::Vector{Int64}`: Multiplicities of unitcells inside
    of a supercell in each of the lattice vector directions. The `j`-th element
    corresponds to the multiplicity in the lattice vector direction given by
    `unit_lattvecs[j,:]`.

# Output

  - `folds::Vector{Int64}`: `unit_point_frac` was the origin of a unitcell in the
    ultracell in fractional coordinates that got folded back into the supercell
    at the origin. `folds` is the result from the folding process whose elements
    sattisfy `1<=folds[j]<=ums[j]` where `ums` is short for `unit_multiplicity_super`.
"""
function fold_to_qe_origin(
    unit_point_frac::Vector{Int64},
    unit_multiplicity_super::Vector{Int64},
)
    # Espresso uses 1-based indexing in their force constants and not 
    # 0-based as in our coordinate system.
    folds = mod.(unit_point_frac .+ 1, unit_multiplicity_super)

    for (i, folding) in enumerate(folds)
        # Smaller-than is just for safety but the mod-function of julia 
        # should be strictly positive definite (not like fortran)
        if folding <= 0
            folds[i] += unit_multiplicity_super[i]
        end
    end

    return folds
end
