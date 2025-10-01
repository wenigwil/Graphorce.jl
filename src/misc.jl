function read_highsympath(file::AbstractString)
    sympathfile = open(file)

    numlines = parse(Int64, readline(sympathfile))

    sympath = Matrix{Float64}(undef, (numlines, 3))
    for line in 1:numlines
        test = parse.(Float64, split(readline(sympathfile)))
        sympath[line, :] = test
    end

    close(sympathfile)
    return sympath
end

"""
    function points_to_path(
        pointlist::Matrix{Float64};
        numpoints_per_section::Int64 = 50,
        return_stitches = false)

Generate a path by parwise-linear interpolation between all points from an array of
points. Control the step size by choosing the number of steps in the linear
interpolation between each point pair (section).

The `return_stitches` argument ensures compatibility with the `path_to_distance()`
function. Every point from `pointlist` will be in the returned path.

# Arguments

  - `pointlist::Matrix{Float64}`: Collection of points which will be interpolated.
    The `j`-th point corresponds to `pointlist[j,:]`.
  - `numpoints_per_section::Int64`: The number of steps from one point to another
    including the arrival point itself.
    `pointlist` in the returned path.
  - `return_stitches::Bool`: Whether or not a `BitVector` should be returned where
    two consecutive `true`-values mark a stichting where multiple `pointlists` are
    joined. (Used for plotting dispersions with a symmetry path containing U|K for
    example)

# Output

  - `path::Matrix{Float64}`: Collection points resulting from the interpolation. The
    size of the resulting array should be
    `( (size(pointlist,1)-1)*numpoints_per_section + 1, 3 )`
  - `stitches::BitVector`: Vector of size of `path` that is `false` everywhere except
    on positions where two `path`s from two `pointlists` are stitched together to form
    one connected path.
"""
function points_to_path(
    pointlist::Matrix{Float64};
    numpoints_per_section::Int64 = 50,
    return_stitches::Bool = false,
)
    numpoints = size(pointlist, 1)
    numsections = numpoints - 1
    numpoints_path = numsections * numpoints_per_section + 1

    path = zeros(Float64, (numpoints_path, 3))
    path[1, :] = pointlist[1, :]
    section = 0
    for is in 1:numsections
        start = 1 + numpoints_per_section * section
        stop = 1 + numpoints_per_section * (section + 1)

        direction = pointlist[is + 1, :] - pointlist[is, :]
        direction ./= numpoints_per_section
        i = 0
        for js in start:stop
            path[js, :] .= i * direction + path[start, :]
            i += 1
        end
        section += 1
    end

    # This part is only for simplifying the path_to_distance function later on
    stitches = BitVector(zeros(Int16, size(path, 1)))
    if return_stitches
        return path, stitches
    else
        return path
    end
end

"""
    function points_to_path(
        points::Matrix{Float64},
        more_pointlists...;
        numpoints_per_section::Int64 = 50,
        return_stitches = false)

Method of the `points_to_path` function that is capable of handling more than one
`pointlist`.
"""
function points_to_path(
    points::Matrix{Float64},
    more_pointlists...;
    numpoints_per_section::Int64 = 50,
    return_stitches = false,
)
    concat_path =
        points_to_path(points; numpoints_per_section = numpoints_per_section)

    stitches = BitVector(zeros(Int16, size(concat_path, 1)))
    for pointlist in more_pointlists
        stitches[end] = 1
        path =
            points_to_path(pointlist, numpoints_per_section = numpoints_per_section)
        new_stitch = BitVector(zeros(Int16, size(path, 1)))
        new_stitch[1] = 1
        concat_path = vcat(concat_path, path)
        stitches = vcat(stitches, new_stitch)
    end

    if return_stitches
        return concat_path, stitches
    else
        return concat_path
    end
end

"""
    function xticks_from_path(
        path::Matrix{Float64},
        sympoints::Matrix{Float64},
        sympoint_labels::Vector{String})

Generate the positions and labels of the ticks on x-axis for plotting dispersions.

# Arguments

  - `path::Matrix{Float64}`: Output from the `points_to_path()`-function.

  - `sympoints::Matrix{Float64}`: A list of symmetry points that is supposedly
    contained in `path`. Must match the ordering of `sympoint_labels`
  - `sympoint_labels::Vector{String}`: A list of the labels for each point in
    `sympoints`. The ordering must match with `sympoints`.

# Output

  - `xtick_labels::Vector{String}`: A collection of the labels for each symmetry
    point found in `path`. Also constructs labels like "U|K" if two symmetry points
    from `sympoints` are found directly next to each other which indicates a
    stitching.

  - `xtick_positions::Vector{Float64}`: A collection of the distances at which an
    xtick-label should be put.
"""
function xticks_from_path(
    path::Matrix{Float64},
    sympoints::Matrix{Float64},
    sympoint_labels::Vector{String},
)
    # Calculate the minimum component step size done in the path calculation
    unique_path = unique(path)
    minstep = minimum(abs.(unique_path[2:end, :] - unique_path[1:(end - 1), :]))
    # Array that marks the sympoint positions in the path array
    sympoint_positions = zeros(Int16, size(path, 1))
    # Array that holds the characters for the xaxis ticks
    for isymp in axes(sympoint_labels, 1)
        sympoint = sympoints[isymp, :]
        for jpath in axes(path, 1)
            if isapprox(path[jpath, :], sympoint; atol = 0.1 * minstep)
                sympoint_positions[jpath] = isymp
            end
        end
    end

    xticks = Vector{String}(undef, 0)
    push!(xticks, sympoint_labels[sympoint_positions[1]])
    for jpath in 2:(size(path, 1) - 1)
        if sympoint_positions[jpath - 1] > 0 && sympoint_positions[jpath] > 0
            label1 = sympoint_labels[sympoint_positions[jpath - 1]]
            label2 = sympoint_labels[sympoint_positions[jpath]]
            label = label1 * "|" * label2
            push!(xticks, label)
        elseif sympoint_positions[jpath] > 0 && sympoint_positions[jpath + 1] == 0
            label = sympoint_labels[sympoint_positions[jpath]]
            push!(xticks, label)
        end
    end
    push!(xticks, sympoint_labels[sympoint_positions[end]])

    return xticks, sympoint_positions
end

"""
    path_to_distance(path::Matrix{Float64}, stitches::BitVector)

Converts a given `path` into a plotable array of accumulated distance. Using the
`stitches` returned by `points_to_path()` the plotable array is generated for paths
that are the product of to paths stitched together.
"""
function path_to_distance(path::Matrix{Float64}, stitches::BitVector)
    numpoints = size(path, 1)

    distances = Vector{Float64}(undef, numpoints)

    distances[1] = 0.0
    for i in 1:(numpoints - 1)
        if stitches[i] + stitches[i + 1] == 2
            distances[i + 1] = distances[i]
        else
            dist = LinAlg.norm(path[i, :] - path[i + 1, :])
            distances[i + 1] = distances[i] + dist
        end
    end

    return distances
end

"""
    path_to_xaxis(
        path::Matrix{Float64},
        stitches::BitVector,
        sympoints::Matrix{Float64},
        sympoint_labels::Vector{String})

Process a symmetry path such that you can plot properties along it. Generate a
plotable array from an array of points and also return the xtick-labels and -positions
of them.
"""
function path_to_xaxis(
    path::Matrix{Float64},
    stitches::BitVector,
    sympoints::Matrix{Float64},
    sympoint_labels::Vector{String},
)
    # First grab the information for the xticks and the positions of the xticks
    xticks_labels, sympoint_positions =
        xticks_from_path(path, sympoints, sympoint_labels)

    # Contract the path to a collection of distances between each point in it
    # Stitches have to be used here to correctly handle things like "U|K"
    distances = path_to_distance(path, stitches)

    # Grab the positions of the xticks at the plotable distances.
    # Make sure the xticks_pos is unique because at the stitching they double up
    xticks_pos = unique(distances[sympoint_positions .> 0])

    return distances, xticks_labels, xticks_pos
end
