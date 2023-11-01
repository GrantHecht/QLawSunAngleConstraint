# Web location of default kernels
ENV["WEB_LSK"] = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
ENV["WEB_SPK"] = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"

# Local locations of default kernels
ENV["LSK"] = joinpath(@__DIR__, "..", "..", "data", "spice", "naif0012.tls")
ENV["SPK"] = joinpath(@__DIR__, "..", "..", "data", "spice", "de440.bsp")

"""
    download_kernals()

Downloads the kernals if they don't exists locally.
"""
function download_kernals()
    isfile(ENV["LSK"]) || download(ENV["WEB_LSK"], ENV["LSK"])
    isfile(ENV["SPK"]) || download(ENV["WEB_SPK"], ENV["SPK"])
    return nothing
end

"""
    furnsh_kernals()

Loads the meta kernal file into the spice system.
"""
function furnsh_kernals()
    download_kernals()
    furnsh(ENV["LSK"])
    furnsh(ENV["SPK"])
    return nothing
end

"""
    unload_kernals()

Unloads the meta kernal file from the spice system.
"""
function unload_kernals()
    unload(ENV["LSK"])
    unload(ENV["SPK"])
    return nothing
end
