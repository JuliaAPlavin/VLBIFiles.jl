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

axes_dicts_all(fh::FITSHeader) = @p axis_types(fh) map(axis_dict(fh, _))

function axis_dict(fh::FITSHeader, ctype::AbstractString)
    ind = axis_ind(fh, ctype)
    re = (@oget fh["XTENSION"]) == "BINTABLE" ?
        Regex("^(MAXIS|CTYPE|CDELT|CRPIX|CRVAL)$(ind)\$") :  # fits idi / bintables in general
        Regex("^([A-Z]+)$(ind)\$")  # uvfits, also hardcode whitelist?
    matching_cards = [match(re, k)[1] => v for (k, v) in pairs(fh) if occursin(re, k)]
    return Dict(matching_cards)
end

axis_val(args...; zero_reference=false) = axis_vals(args...; zero_reference) |> only
axis_vals(fh::FITSHeader, ctype; zero_reference=false) = axis_vals(axis_dict(fh, ctype); zero_reference)
function axis_vals(dict::Dict; zero_reference=false)
    n = @oget dict["NAXIS"] dict["MAXIS"]
    if zero_reference
        return ((1:n) .- dict["CRPIX"]) .* dict["CDELT"]
    else
        return dict["CRVAL"] .+ ((1:n) .- dict["CRPIX"]) .* dict["CDELT"]
    end
end


function named_axiskeys_tablecol(fh::FITSHeader)
    ax_dicts = axes_dicts_all(fh)
    @p ax_dicts map() do r
        name = Symbol(r["CTYPE"])
        name => axis_vals(r)
    end NamedTuple
end

const val_to_stokes = Dict(-4 => :LR, -3 => :RL, -2 => :LL, -1 => :RR, 1 => :I, 2 => :Q, 3 => :U, 4 => :V)

function named_axiskeys_tablecol_fitsidi(fh::FITSHeader)
    base = named_axiskeys_tablecol(fh)
    valmaps = (COMPLEX=(@o [:re, :im, :wt][Int(_)]), STOKES=(@o val_to_stokes[_]))
    updvals = map(keys(base), values(base)) do k, origval
        valmap = @oget valmaps[k]
        isnothing(valmap) && return origval
        return valmap.(origval)
    end
    NamedTuple{keys(base)}(updvals)
end
