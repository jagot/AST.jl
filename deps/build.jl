paths_file = "$(Pkg.dir("AST"))/src/paths.jl"

if !isfile(paths_file)
    atsp,grasp = try
        ENV["ATSP"], ENV["GRASP"]
    catch
        error("Did not find ATSP and/or GRASP environment variables. Please define them and rebuild AST.")
    end

    isfile("$atsp/hf") || error("Could not find hf exe in $atsp")
    isfile("$grasp/rmcdhf") || error("Could not find rmcdhf exe in $grasp")

    println("Found ATSP2K at $atsp and GRASP2K at $grasp")

    open(paths_file, "w") do file
        write(file, """atsp = \"$atsp\"
    grasp = \"$grasp\"
    """)
    end
end
