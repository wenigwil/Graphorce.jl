#=
  ____                                    
 |  _ \  __ _  _ __  ___   ___  _ __  ___
 | |_) |/ _` || '__|/ __| / _ \| '__|/ __|
 |  __/| (_| || |   \__ \|  __/| |   \__ \
 |_|    \__,_||_|   |___/ \___||_|   |___/                        
==========================================

This file contains the Parsers module in which objects for holding
data are defined in! Each object will correspond to a file from an
external application such as elphbolt, Quantum Espresso or similar.
The constructors (inner and outer) of the objects are used to parse
the file that they belong to. This is mostly for convenience when
accessing data and modularization...
=#

module Parsers

import Logging

struct ebInputData
    input_file::AbstractString
    allocations::AbstractDict
    crystal_info::AbstractDict
    numerics::AbstractDict

    function ebInputData(
        # Input file that will be read in
        input_file::AbstractString,
        # elphbolt namelists contain quantities that need to be initialized 
        # with dimensions that depend on the numelements and numatoms
        numelements::Integer,
        numatoms::Integer,
        # TODO: For know this struct implements the allocations and 
        # crystal_info namelists from elphbolt. crystal_info depends in 
        # at least 6 quantities on numelements or numatoms.
        # The electrons namelist depends on numconc and numT! Implement 
        # those
    )
        allocations = Dict{String,Any}()
        crystal_info = Dict{String,Any}()
        numerics = Dict{String,Any}()

        # Initialize the allocations namelist with the numelements and 
        # numatoms that have to be passed in any way
        Logging.@info "Inner ebID-Constructor: Initializing empty allocations dictionary..."
        init_allocations!(allocations, numelements, numatoms)

        # Initialize the crystal_info namelist. Quantities in this namelist 
        # depend strictly on numelements and numatoms
        Logging.@info "Inner ebID-Constructor: Initializing empty crystal_info dictionary..."
        init_crystal_info!(crystal_info, numelements, numatoms)

        # Initialize the numerics namelist. This namelist does not depend 
        # on numelements or numatoms
        Logging.@info "Inner ebID-Constructor: Initializing empty numerics dictionary..."
        init_numerics!(numerics)

        # Construct the ebInputData-Object
        new(input_file, allocations, crystal_info, numerics)
    end
end

function ebInputData(input_file_path::AbstractString)
    # Check that the input_file actually exists
    if isfile(input_file_path) == false
        Logging.@error "Outer ebID-Constructor: Given path is not a regular file!"
        return
    else
        Logging.@debug "Outer ebID-Constructor: Given path is a regular file" input_file_path
    end

    # First read the whole file into a vector of strings
    file_vector = readlines(input_file_path)
    Logging.@debug "Outer ebID-Constructor: Read the file into a vector of strings..."

    # Cleaning of blank lines and lines of beginning /
    clean_mask = BitMatrix(undef, file_vector.size[1], 1)
    for (i, line) in enumerate(file_vector)
        # Matching regex expressions against each line
        if match(r"^/", line) !== nothing
            # Mark the mask for lines beginning with /
            clean_mask[i] = true
        elseif match(r"^$", line) !== nothing
            #Mark the mask for blank lines
            clean_mask[i] = true
        end
    end
    # Invert the mask
    clean_mask = .!clean_mask

    # Strip file_vector from the lines marked in clean_mask
    file_vector = file_vector[clean_mask]
    Logging.@debug "Outer ebID-Constructor: Cleaned empty lines and lines beginning with \"/\"..."

    # Create a mask that marks namelist headers (allocations, crystal_info 
    # ...) and indicates to which namelist each attribute belongs
    cut_mask = Vector{String}(undef, file_vector.size[1])
    nml_header = String
    for (i, line) in enumerate(file_vector)
        # Match the namelist header with regex and otherwise mark the mask 
        # with the header name (nml_name)
        if match(r"^&", line) !== nothing
            # Double split for whitespace-safety
            nml_header = split(split(line, "&")[2])[1]
            # Differentiate the header from the attributes
            cut_mask[i] = "head " * nml_header
        else
            cut_mask[i] = nml_header
        end
    end

    Logging.@debug "Outer ebID-Constructor: Created an address mask for each namelist..." cut_mask
    # Initiate and set the two attributes that all namelists depend on so 
    # we can call the inner constructor with these
    # This is basically going through the allocations namelist and put it 
    # into the corresponding dictionary in the ebInputData-object
    numelements = 0
    numatoms = 0
    # Loop through the part of the file that corresponds to the allocations
    Logging.@info "Outer ebID-Constructor: Starting to read the allocations-namelist to create an empty ebInputData-object... "
    for line in file_vector[cut_mask .== "allocations"]
        nml_attr = split(split(line, "=")[1])[1]
        nml_val = parse(Int64, split(split(line, "=")[2])[1])
        # Because the inner constructor is not called yet we have no access 
        # to the types of each nml_attribute yet and have to do this crude 
        # checking
        if nml_attr == "numelements"
            numelements = nml_val
        elseif nml_attr == "numatoms"
            numatoms = nml_val
        else
            # This will be thrown if the format or names in the allocations 
            # namelist change!
            Logging.@error "Outer ebID-Constructor: Attribute defined in allocations namelist that is not valid!"
        end
    end

    # We have read the numelements and numatoms so we can now initialize an 
    # ebInputData-object with the default values and array sizes that 
    # depend on numelements and numatoms
    dataobject = ebInputData(input_file_path, numelements, numatoms)
    Logging.@debug "Outer ebID-Constructor: Empty ebInputData-object initialized with " numelements numatoms input_file_path

    read_namelist_into_dict!(file_vector, cut_mask, "crystal_info", dataobject.crystal_info)
    read_namelist_into_dict!(file_vector, cut_mask, "numerics", dataobject.numerics)

    return dataobject
end

function read_namelist_into_dict!(
    file_vector::Vector{String},
    addr_mask::Vector{String},
    nml_header::String,
    dict::AbstractDict,
)
    Logging.@info "Outer ebID-Constructor: Starting to read the $(nml_header)-namelist into an ebID-object... "
    for line in file_vector[addr_mask .== nml_header]
        # Crudely disect the line into attribute, (index, )value
        disect = disect_nml_line(line)
        attr = disect[1]

        # Get the type from the attribute out of disect
        ele_type, conc_type = get_attr_type(dict, attr)

        if isa(dict[attr], Array)
            # Convert the array of string ranges into an 
            # Array{UnitRange{Int64}}
            ranges = rangestrings_to_unitranges(disect[3], dict[attr])

            clean_vals = clean_array_val(dict[attr], disect[2])

            try
                dict[attr][ranges...] = parse.(ele_type, clean_vals)
            catch
                dict[attr][ranges...] = convert.(ele_type, clean_vals)
            end
            Logging.@debug """
            Outer ebID-Constructor: Disected and cleaned a line 
            """ line attr conc_type ranges clean_vals
        else
            clean_val = clean_non_array_val(dict[attr], disect[2])

            # Try to parse the value string into the concrete type. If that 
            # does not work, the corresponding attribute is of type string 
            # and then it has to be converted
            try
                dict[attr] = parse(conc_type, clean_val)
            catch
                dict[attr] = convert(conc_type, clean_val)
            end
            Logging.@debug """
            Outer ebID-Constructor: Disected and cleaned a line 
            """ line attr ele_type conc_type clean_val
        end
    end
end

function rangestrings_to_unitranges(
    indices::Vector{SubString{String}},
    indexed_array::Array,
)
    ranges = Array{UnitRange{Int64}}(undef, length(size(indexed_array)))

    # Special case for 1d attributes that dont need indices to be specified
    if size(indices) == (0,) && length(ranges) == 1
        ranges[1] = 1:size(indexed_array)[1]
    end

    for (i, range_string) in enumerate(indices)
        # A colon specifies the whole range of the indexed array in that 
        # dimension
        if range_string == ":"
            ranges[i] = 1:size(indexed_array)[i]
        else
            ranges[i] = parse(Int64, range_string):parse(Int64, range_string)
        end
    end
    return ranges
end

function clean_array_val(dict_attr::Any, array_val::AbstractString)
    values = split(array_val)
    # Clean the array
    for (i, value) in enumerate(values)
        values[i] = clean_non_array_val(dict_attr, value)
    end
    return values
end

# TODO: does dict_attr have to be any?
function clean_non_array_val(dict_attr::Any, val::AbstractString)
    clean_val = val
    if isa(dict_attr, String) || isa(dict_attr, AbstractArray{String})
        # \"Si\" --> Si
        clean_val = strip(val, ['"', ''', '.', ' '])
    elseif isa(dict_attr, Bool) || isa(dict_attr, AbstractArray{Bool})
        # .true. --> true
        clean_val = strip(val, ['.', ' '])
    end

    return clean_val
end

function get_attr_type(nml_dict::AbstractDict, attr::AbstractString)
    # Full type (as in Array{Float64})
    conc_type = typeof(nml_dict[attr])
    # Type of the elements (in the example above Float64)
    ele_type = eltype(nml_dict[attr])
    return ele_type, conc_type
end

function disect_nml_line(line::AbstractString)
    # First split away the comment
    noco = split(line, "!")[1]

    # Split name-side and value-side
    name = split(split(noco, "=")[1])[1]
    val = split(noco, "=")[2]

    # Split the name into attribute and index
    attr_and_index = split(name, "(")
    attr = attr_and_index[1]
    if attr_and_index.size[1] > 1
        range = chop(split(name, "(")[2], tail = 1)

        #Split index into indices
        indices = split(range, ",")
        return [attr, val, indices]
    else
        return [attr, val, Array{SubString{String}}(undef, 0)]
    end
end

function init_allocations!(
    allocations::AbstractDict,
    numelements::Integer,
    numatoms::Integer,
)
    allocations["numelements"] = numelements
    allocations["numatoms"] = numatoms
end

function init_crystal_info!(ci::AbstractDict, numelements::Integer, numatoms::Integer)
    ci["name"] = "Crystal"
    ci["elements"] = Vector{String}(undef, numelements)
    ci["atomtypes"] = Vector{Int64}(undef, numatoms)
    ci["masses"] = Vector{Float64}(undef, numelements)
    ci["VCA"] = false
    ci["DIB"] = false
    ci["lattvecs"] = Array{Float64}(undef, (3, 3))
    ci["basis"] = Array{Float64}(undef, (3, numatoms))
    ci["polar"] = false
    ci["born"] = Array{Float64}(undef, (3, 3, numatoms))
    ci["epsilon"] = Array{Float64}(undef, (3, 3))

    # This one is labeled on the elphbolt github page as a Real
    ci["read_epsiloninf"] = false

    ci["epsiloninf"] = 0.0
    ci["epsilon0"] = 0.0
    ci["T"] = -1.0
    ci["twod"] = false
    ci["subs_masses"] = Vector{Float64}(undef, numelements)
    ci["subs_conc"] = Vector{Float64}(undef, numelements)
    ci["bound_length"] = 1e12
    ci["specfac"] = 0.0
end

function init_numerics!(nu::AbstractDict)
    nu["qmesh"] = Vector{Int64}(undef, 3)
    nu["mesh_ref"] = 1
    nu["fsthick"] = 0.0
    nu["datadumpdir"] = "./"
    nu["read_gq2"] = false
    nu["read_gk2"] = false
    nu["read_V"] = false
    nu["W_OTF"] = true
    nu["tetrahedra"] = false
    nu["fourph"] = false
    nu["four_mesh_ref"] = 1
    nu["phe"] = false
    nu["Y_OTF"] = true
    nu["phise"] = false
    nu["phiso_1b_theory"] = "DIB-1B"
    nu["phsubs"] = false
    nu["phbound"] = false
    nu["onlyphbte"] = false
    nu["elchimp"] = false
    nu["elbound"] = false
    nu["onlyebte"] = false
    nu["drag"] = true
    nu["maxiter"] = 50
    nu["conv_thres"] = 1e-4
    nu["runlevel"] = 1
    nu["plot_along_path"] = false
    nu["ph_en_min"] = 0.0
    nu["ph_en_max"] = 1.0
    nu["el_en_min"] = 0.0
    nu["el_en_max"] = 1.0
    nu["el_en_num"] = 100
    nu["ph_mfp_npts"] = 100
    nu["ph_abs_q_npts"] = 100
    nu["use_Wannier_ifc2s"] = false
    nu["solve_bulk"] = true
    nu["solve_nano"] = false
end

"""
    dfpt_qeOutputData(output_file_path::String)

Type for the output data of a Quantum Espresso DFPT
calculation that is read-in from an output file at output_file_path.
"""
struct dfpt_qeOutputData
    output_file::AbstractString
    properties::Dict

    function dfpt_qeOutputData(output_file_path::String)
        # Check that file exists
        if isfile(output_file_path) == false
            Logging.@error "Inner dfpt_qeOD-constructor: Given path is not a regular file!"
            return
        else
            Logging.@debug "Inner dfpt_qeOD-constructor: Given path is a regular file" output_file_path
        end

        # Reading the header
        output_file = open(output_file_path)

        # FIRST LINE
        # Reading number of species (ntype), number of atoms (nat), type of 
        # bravais lattice (?) (ibrav), and a vector of cell dimensions (celldm)
        celldm = Vector{Float64}(undef, 6)
        ntype,
        nat,
        ibrav,
        celldm[1],
        celldm[2],
        celldm[3],
        celldm[4],
        celldm[5],
        celldm[6] = parse.(Float64, split(readline(output_file)))
        # Convert to Ints so we can loop over them
        ntype = convert(Int64, ntype)
        nat = convert(Int64, nat)
        ibrav = convert(Int64, ibrav)

        Logging.@debug "readifc2() First line read..." nat ntype ibrav celldm

        # TODO Deal with the special ibrav=0 case
        # For this we have to find out what quantum espresso prints into the 
        # ifc2-output file in that case... 15 year old forums may suggest a 
        # solution lol

        # NEXT ntype LINE(S)
        # ntype are the number species of atoms
        mass = Vector{Float64}(undef, ntype)
        label = Array{String}(undef, ntype)
        for i in 1:ntype
            temp = split(readline(output_file))
            label[i] = chop(temp[2], head = 1, tail = 1)
            mass[i] = parse(Float64, temp[3])
        end

        # Convert the masses from the Rydberg-ian unitsystem used in 
        # Quantum-Espresso into *relative atomic mass*

        Logging.@debug "readifc2() Next $ntype line(s) read..." mass label

        # NEXT nat LINE(S)
        # From these lines we get which species (tipo) sits on which basis 
        # position (basis_pos)
        tipo = Vector{Int64}(undef, nat)
        basis_pos = Matrix{Float64}(undef, (nat, 3))
        for i in 1:nat
            temp = split(readline(output_file))
            tipo[i] = parse(Int64, temp[2])
            basis_pos[i, :] = parse.(Float64, temp[3:end])
        end

        Logging.@debug "readifc2() Next $nat line(s) read..." tipo basis_pos

        # NEXT LINE
        # Is the material a polar one?
        polar_key = split(readline(output_file))[1]

        Logging.@debug "readifc2() Next 1 line read..." polar_key

        if polar_key == "T"
            # NEXT 3 LINES
            # Read the dielectric tensor
            dielectric_tensor = Matrix{Float64}(undef, (3, 3))
            for i in 1:3
                dielectric_tensor[i, :] = parse.(Float64, split(readline(output_file)))
            end
            Logging.@debug "readifc2() Next 3 lines read... Dielectric Tensor" dielectric_tensor

            # NEXT 4*nat LINES
            # Read the Born effective charges
            zeff = Array{Float64}(undef, (nat, 3, 3))
            for i in 1:nat
                # We are not reading the first index. We already know that
                readline(output_file)
                for j in 1:3
                    zeff[i, j, :] = parse.(Float64, split(readline(output_file)))
                end
            end
            Logging.@debug "readifc2() Next $(4*nat) lines read... Born-Eff. Charges" zeff
        end

        # NEXT LINE
        # Read the super cell sizes
        qcell = parse.(Int32, split(readline(output_file)))
        Logging.@debug "readifc2() Next 1 line read... Super Cell Size" qcell

        # NEXT 3*3*nat*nat*( qcell[1]*qcell[2]*qcell[3] + 1 ) lines
        # First we define the number of force-constant-supercell-fields 
        nfc2 = 3 * 3 * nat * nat
        no_qcells = qcell[1] * qcell[2] * qcell[3]

        # Allocating the force tensor
        ifc2_tensor_dims = (3, 3, nat, nat, qcell[1], qcell[2], qcell[3])
        ifc2_tensor = Array{Float64}(undef, ifc2_tensor_dims)

        Logging.@debug """
        Allocated the force-constant tensor with following dimensions
        ifc2_tensor_dims = (3, 3, nat, nat, qcell[1], qcell[2], qcell[3])\n
        """ ifc2_tensor_dims

        # Reading the force constant-tensor.
        for i in 1:nfc2
            # Reading the atom displacements and atom/species addresses
            i_displ, j_displ, i_at_addr, j_at_addr =
                parse.(Int64, split(readline(output_file)))

            # Reading the supercell coordinates and actual ifc2-tensor elements
            for j in 1:no_qcells
                temp = split(readline(output_file))

                i_qcell, j_qcell, k_qcell = parse.(Int64, temp[1:(end - 1)])
                ifc2_tensor[
                    i_displ,
                    j_displ,
                    i_at_addr,
                    j_at_addr,
                    i_qcell,
                    j_qcell,
                    k_qcell,
                ] = parse(Float64, temp[end])
            end
        end
        close(output_file)

        # Create the properties dictionary for storing the read-in data
        properties = Dict{String,Any}()

        # Assign the read-in data to the properties dictionary
        properties["numspecies"] = ntype
        properties["numbasisatoms"] = nat
        properties["ibrav"] = ibrav
        properties["celldm"] = celldm
        properties["species2masses"] = mass
        properties["species"] = label
        properties["basisatom2species"] = tipo
        properties["basis"] = basis_pos
        properties["polar"] = polar_key
        properties["epsilon"] = dielectric_tensor
        properties["born"] = zeff
        properties["qcell"] = qcell
        properties["ifc2"] = ifc2_tensor

        # Advise the constructor to return a dfpt_qeOutputData-object
        new(output_file_path, properties)
    end
end

# End of module Parsers
end
