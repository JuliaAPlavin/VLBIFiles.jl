module RectiGridsAxisKeysExtraExt

using RectiGrids
using AxisKeysExtra
using VLBIFiles
using VLBIFiles: FitsImage, SVector

function VLBIFiles.image_clean(::Type{KeyedArray}, fimg::FitsImage; kwargs...)
    @eval InterferometricModels.exp_for_gaussintensity(x) = InterferometricModels.exp_for_gaussintensity_difmap(x)
    result = try
        @invokelatest _image_clean(KeyedArray, fimg; kwargs...)
    finally
        @eval InterferometricModels.exp_for_gaussintensity(x) = InterferometricModels.exp_for_gaussintensity_basic(x)
    end
    return result
end

function _image_clean(::Type{KeyedArray}, fimg::FitsImage; kwargs...)
    im = VLBI.image_clean(fimg; kwargs...)
    return im.(with_axiskeys(grid)(SVector, fimg.data))
end

end
