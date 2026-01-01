module RectiGridsAxisKeysExtraExt

using RectiGrids
using AxisKeysExtra
using VLBIFiles
using VLBIFiles: FitsImage, SVector

function VLBIFiles.image_clean(::Type{KeyedArray}, fimg::FitsImage; kwargs...)
    im = VLBI.image_clean(fimg; kwargs...)
    return im.(with_axiskeys(grid)(SVector, fimg.data))
end

end
