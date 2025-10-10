
function read_ifc3(output_file_path::String)
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

    close(output_file)
end
