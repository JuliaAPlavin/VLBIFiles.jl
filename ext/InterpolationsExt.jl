module InterpolationsExt

using Interpolations
using VLBIFiles

VLBIFiles.image_stored(fimg::VLBI.FitsImage) = linear_interpolation(VLBI.image_stored(VLBIFiles.KeyedArray, fimg); extrapolation_bc=zero(eltype(fimg.data)))

end
