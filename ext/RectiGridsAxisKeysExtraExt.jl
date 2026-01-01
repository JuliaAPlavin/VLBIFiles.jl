module RectiGridsAxisKeysExtraExt

using RectiGrids
using AxisKeysExtra
using VLBIFiles
using VLBIFiles: FitsImage, SVector

function VLBIFiles.image_clean(::Type{KeyedArray}, fimg::FitsImage; kwargs...)
    @eval InterferometricModels.exp_for_gaussintensity(x) = InterferometricModels.exp_for_gaussintensity_difmap(x)
    result = try
        im = @invokelatest VLBI.image_clean(fimg; kwargs...)
        im.(with_axiskeys(grid)(SVector, fimg.data))
    finally
        @eval InterferometricModels.exp_for_gaussintensity(x) = InterferometricModels.exp_for_gaussintensity_basic(x)
    end
    return result
end

end
