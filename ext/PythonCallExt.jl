module PythonCallExt

using PythonCall
using VLBIFiles: VLBIFiles, UVData

function VLBIFiles.read_data_raw(uvdata::UVData, ::typeof(pyimport))
    np = pyimport("numpy")
    astropy_fits = pyimport("astropy.io.fits")

    builtins = pyimport("builtins")
    hdus = astropy_fits.open(builtins.open(uvdata.path, "rb"), memmap=false)
    raw = hdus[0].data
    hdus.close()

    names = pyconvert(Vector{String}, raw.dtype.names)
    result = Dict{Symbol,Any}()
    for n in names
        col = raw[n]
        native = np.ascontiguousarray(col.astype(col.dtype.newbyteorder("=")))
        result[Symbol(n)] = pyconvert(Array, native)
    end
    result[:DATA] = permutedims(result[:DATA], reverse(1:ndims(result[:DATA])))
    return result
end

end
