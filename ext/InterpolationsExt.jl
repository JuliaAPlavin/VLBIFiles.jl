module InterpolationsExt

using Interpolations
using VLBIFiles

VLBIFiles.image_stored(fimg::VLBI.FitsImage) = linear_interpolation(VLBI.image_stored(VLBIFiles.KeyedArray, fimg); extrapolation_bc=zero(eltype(fimg.data)))

# better to subtract and then interpolate, than the other way around
VLBIFiles.image_residual(fimg::VLBI.FitsImage) = linear_interpolation(VLBI.image_residual(VLBIFiles.KeyedArray, fimg); extrapolation_bc=zero(eltype(fimg.data)))

end
