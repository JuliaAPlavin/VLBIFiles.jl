function lazycolumntable(hdu::TableHDU)
    _ensure_hdu_active!(hdu)
    colnames = columnnames(hdu) |> Tuple
    ctx = _mmap_table_context(hdu)

    cols = if ctx !== nothing
        map(colnames) do colname
            colnum = fits_get_colnum(hdu.fitsfile, String(colname); case_sensitive=true)
            col = _mmap_column(ctx, hdu, colnum)
            col !== nothing ? col : TableHDUColumn(hdu, colname)
        end
    else
        map(colname -> TableHDUColumn(hdu, colname), colnames)
    end
    tbl = NamedTuple{colnames}(cols)

    fhead = FITSIO.read_header(hdu)
    haskey(tbl, :FLUX) && @reset tbl.FLUX = mapview(colarray_postfunc(fhead, Val(:FLUX)), tbl.FLUX)
    haskey(tbl, :WEIGHT) && @reset tbl.WEIGHT = mapview(colarray_postfunc(fhead, Val(:WEIGHT)), tbl.WEIGHT)

    return tbl
end

function colarray_postfunc(fhead::FITSHeader, ::Val{:FLUX})
    col_naxkeys = named_axiskeys_tablecol_fitsidi(fhead)
    x -> KeyedArray(reshape(x, Tuple(map(length, col_naxkeys))); col_naxkeys...)
end

function colarray_postfunc(fhead::FITSHeader, ::Val{:WEIGHT})
    col_naxkeys = named_axiskeys_tablecol_fitsidi(fhead)[(:STOKES, :BAND)]
    resf(x::AbstractArray) = KeyedArray(reshape(x, Tuple(map(length, col_naxkeys))); col_naxkeys...)
	resf(x::Number) = x
	return resf
end

function colarray_postfunc(fhead::FITSHeader, ::Val{:DATA})
    col_naxkeys = named_axiskeys_tablecol_uvfits(fhead)
    x -> KeyedArray(reshape(x, Tuple(map(length, col_naxkeys))); col_naxkeys...)
end



struct TableHDUColumn{T, P} <: AbstractVector{T}
    hdu::TableHDU
    properties::P
end

function TableHDUColumn(hdu::TableHDU, colname::Union{AbstractString,Symbol})
    colname = String(colname)
    _ensure_hdu_active!(hdu)
    nrows = fits_get_num_rows(hdu.fitsfile)
    colnum = fits_get_colnum(hdu.fitsfile, colname; case_sensitive=true)
    T, rowsize, isvariable = fits_get_col_info(hdu.fitsfile, colnum)
    ET = isempty(rowsize) ? T : Array{T, length(rowsize)}
    properties = (;name=colname, colnum, length=nrows, rowsize=Tuple(rowsize), isvariable)
    TableHDUColumn{ET, typeof(properties)}(hdu, properties)
end

function _ensure_hdu_active!(hdu::TableHDU)
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
end

Base.size(col::TableHDUColumn) = (col.properties.length,)

function Base.getindex(col::TableHDUColumn{<:Union{Number,AbstractString}}, i::Int)
    _ensure_hdu_active!(col.hdu)
    A = Vector{eltype(col)}(undef, 1)
    CFITSIO.fits_read_col(col.hdu.fitsfile, col.properties.colnum, i, 1, A)
    return only(A)
end

function Base.getindex(col::TableHDUColumn{<:AbstractArray}, i::Int)
    _ensure_hdu_active!(col.hdu)
    A = Vector{eltype(eltype(col))}(undef, col.properties.rowsize)
    CFITSIO.fits_read_col(col.hdu.fitsfile, col.properties.colnum, i, 1, A)
    return A
end

function Base.getindex(col::TableHDUColumn{<:Union{Number,AbstractString}}, I::UnitRange{Int})
    _ensure_hdu_active!(col.hdu)
    A = Vector{eltype(col)}(undef, length(I))
    CFITSIO.fits_read_col(col.hdu.fitsfile, col.properties.colnum, first(I), 1, A)
    return A
end

Base.copyto!(dest::Vector, col::TableHDUColumn{<:Union{Number,AbstractString}}) = CFITSIO.fits_read_col(col.hdu.fitsfile, col.properties.colnum, 1, 1, dest)
