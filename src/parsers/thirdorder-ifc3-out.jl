struct Ifc3Output
    output_file::AbstractString
    properties::Dict{String,Any}

    function Ifc3Output(output_file_path::String)
        if isfile(output_file_path) == false
            Logging.@error """
            read_ifc3: Given path is not a regular file!
            """ output_file_path
            return
        else
            Logging.@debug """
            read_ifc3: Given path is a regular file output_file_path
            """
        end

        output_file = open(output_file_path)

        # Number of 3-atom combinations (triplets)
        # This number of combinations was reduced by symmetry operations
        num_triplets = parse(Int64, readline(output_file))

        # ifc3-Tensor
        ifc3_tensor = Array{Float64}(undef, (num_triplets, 3, 3, 3))

        # Positions of the second and third atoms that are not put into the origin
        trip2position_j = Matrix{Float64}(undef, (num_triplets, 3))
        trip2position_k = Matrix{Float64}(undef, (num_triplets, 3))

        # Atomic indices
        trip2atomindices = Matrix{Int64}(undef, (num_triplets, 3))

        # Skip one line
        readline(output_file)

        for itrip in 1:num_triplets
            # skip the triplet index
            readline(output_file)

            # Positions
            trip2position_j[itrip, :] = parse.(Float64, split(readline(output_file)))
            trip2position_k[itrip, :] = parse.(Float64, split(readline(output_file)))

            # Atomic indices
            trip2atomindices[itrip, :] = parse.(Int64, split(readline(output_file)))

            # Displacement indices and value 
            for idisp in 1:3
                for jdisp in 1:3
                    for kdisp in 1:3
                        tmp = split(readline(output_file))
                        ifc3_tensor[itrip, idisp, jdisp, kdisp] =
                            parse(Float64, tmp[4])
                    end
                end
            end

            # Blank line separation
            readline(output_file)
        end
        close(output_file)

        props = Dict{String,Any}()

        props["ifc3_tensor"] = ifc3_tensor
        props["numtriplets"] = num_triplets
        props["trip2position_j"] = trip2position_j
        props["trip2position_k"] = trip2position_k
        props["trip2atomindices"] = trip2atomindices

        new(output_file_path, props)
    end
end
