using FITSIO: fits_get_num_rows, TableHDU, assert_open, fits_movabs_hdu, fits_get_colnum, fits_get_col_info, CFITSIO, fits_try_read_extname

function name(hdu::FITSIO.HDU)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    return fits_try_read_extname(hdu.fitsfile)
end

axis_types(fh::FITSHeader) = @p begin
    pairs(fh)
    filtermap() do (k, v)
        occursin(r"^CTYPE\d+$", k) ? v : nothing
    end
end

_ismatch(needle::AbstractString, haystack::AbstractString) = needle == haystack
_ismatch(needle::Regex, haystack::AbstractString) = occursin(needle, haystack)
_ismatch(needle, haystack) = false

function axis_ind(fh::FITSHeader, ctype)
    matching_card = @p pairs(fh) filteronly() do (k, v)
        _ismatch(ctype, v)
    end
    key = matching_card.first
    m = match(r"^CTYPE(\d)$", key)
    return parse(Int, m[1])
end

function axis_dict(fh::FITSHeader, ctype::AbstractString)
    ind = axis_ind(fh, ctype)
    re = fh["XTENSION"] == "BINTABLE" ?
        Regex("^(MAXIS|CTYPE|CDELT|CRPIX|CRVAL)$(ind)\$") :  # fits idi / bintables in general
        Regex("^([A-Z]+)$(ind)\$")  # uvfits, also hardcode whitelist?
    matching_cards = [match(re, k)[1] => v for (k, v) in pairs(fh) if occursin(re, k)]
    return Dict(matching_cards)
end

function axis_vals(fh::FITSHeader, ctype; zero_reference=false)
    dict = axis_dict(fh, ctype)
    n = @oget dict["NAXIS"] dict["MAXIS"]
    if zero_reference
        return ((1:n) .- dict["CRPIX"]) .* dict["CDELT"]
    else
        return dict["CRVAL"] .+ ((1:n) .- dict["CRPIX"]) .* dict["CDELT"]
    end
end
    
function axis_val(fh::FITSHeader, ctype; zero_reference=false)
    vals = axis_vals(fh, ctype, zero_reference=zero_reference)
    @assert length(vals) == 1
    return first(vals)
end
