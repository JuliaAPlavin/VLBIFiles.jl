function _geodetic_to_ecef(lat, lon, elev)
    a = 6378137.0        # WGS84 semi-major axis (m)
    f = 1/298.257223563  # WGS84 flattening
    e2 = 2f - f^2
    sinlat, coslat = sincos(lat)
    sinlon, coslon = sincos(lon)
    N = a / sqrt(1 - e2 * sinlat^2)
    x = (N + elev) * coslat * coslon
    y = (N + elev) * coslat * sinlon
    z = (N * (1 - e2) + elev) * sinlat
    SVector(x, y, z)
end

function _parse_mount_type(s)
    ismissing(s) && return VLBIData.AntennaMountType.Unknown
    s == "ALT-AZ" && return VLBIData.AntennaMountType.AltAzimuth
    s == "ALT-AZ+NASMYTH-R" && return VLBIData.AntennaMountType.NaismithR
    s == "ALT-AZ+NASMYTH-L" && return VLBIData.AntennaMountType.NaismithL
    return VLBIData.AntennaMountType.Unknown
end

function _parse_poltypes(s)
    ismissing(s) && return (:U, :U)
    s == "circular" && return (:R, :L)
    s == "linear" && return (:X, :Y)
    return (:U, :U)
end

function ngehtsim_antenna_catalog(path=joinpath(@__DIR__, "data", "Telescope_Site_Matrix.csv"))
    rows = CSV.File(path)
    map(rows) do row
        Antenna(;
            name=Symbol(row.Name),
            xyz=_geodetic_to_ecef(deg2rad(row.Latitude), deg2rad(row.Longitude), row.Elevation),
            mount_type=_parse_mount_type(row.Mount_type),
            poltypes=_parse_poltypes(row.Polarization_basis),
        )
    end
end
