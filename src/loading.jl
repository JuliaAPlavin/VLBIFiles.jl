function guess_type(src)
    src = abspath(src)  # for RFC.File
    line = first(eachline(src))
    if isvalid(line) && startswith(line, r"! \w+")
        MultiComponentModel
    elseif isvalid(line) && startswith(line, "* This file processed by ")
        Alist
    else
        hdunames = map(name, FITS(src))
        ctypes = try
            FITS(src) do f
                read_header(f[1]) |> axis_types
            end
        catch e
            error("Cannot read $src as a FITS file")
        end
        if "RA---SIN" ∈ ctypes
            return FitsImage
        elseif "COMPLEX" ∈ ctypes
            return UVData  # UVFITS
        elseif "UV_DATA" ∈ hdunames
            return UVData  # FITSIDI
        else
            error("Cannot guess data type using FITS header of $src")
        end
    end
end

load(src; kwargs...) = load(guess_type(src), src; kwargs...)
