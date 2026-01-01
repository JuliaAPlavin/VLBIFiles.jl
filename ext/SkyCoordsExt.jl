module SkyCoordsExt

using VLBIFiles
using SkyCoords
using Unitful, UnitfulAngles

"""
    ICRSCoords(uvheader::VLBIFiles.UVHeader)

Extract coordinates of the observed source from a UVHeader object.
"""
function SkyCoords.ICRSCoords(uvheader::VLBIFiles.UVHeader)
    fits_header = uvheader.fits
    ra_rad = deg2rad(fits_header["OBSRA"])
    dec_rad = deg2rad(fits_header["OBSDEC"])
    return ICRSCoords(ra_rad, dec_rad)
end

"""
    ICRSCoords(uvdata::VLBIFiles.UVData)

Extract coordinates of the observed source from a UVData object.
"""
SkyCoords.ICRSCoords(uvdata::VLBIFiles.UVData) = ICRSCoords(uvdata.header)

"""
    ICRSCoords(fitsimage::VLBIFiles.FitsImage)

Extract coordinates of the observed source from a FitsImage object.
"""
function SkyCoords.ICRSCoords(fitsimage::VLBIFiles.FitsImage)
    fits_header = fitsimage.header
    ra_rad = deg2rad(fits_header["OBSRA"])
    dec_rad = deg2rad(fits_header["OBSDEC"])
    return ICRSCoords(ra_rad, dec_rad)
end

end
