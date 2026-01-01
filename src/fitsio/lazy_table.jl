function lazycolumntable(hdu::TableHDU)
	colnames = columnnames(hdu) |> Tuple
	cols = map(colnames) do colname
		TableHDUColumn(hdu, colname)
	end
	NamedTuple{colnames}(cols)
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

# not certain whether batch reads are actually faster...
ALLOW_SCALAR = Ref(true)
allow_scalar(enable::Bool=true) = (ALLOW_SCALAR[] = enable)
function allow_scalar(func::Function, enable::Bool=true)
    old = ALLOW_SCALAR[]
    try
        allow_scalar(enable)
        func()
    finally
        allow_scalar(old)
    end
end
        

function Base.getindex(col::TableHDUColumn{<:Union{Number,AbstractString}}, i::Int)
    ALLOW_SCALAR[] || error("Scalar access to FITS columns is disabled. Enable it by calling `allow_scalar()`.")
    _ensure_hdu_active!(col.hdu)
    A = Vector{eltype(col)}(undef, 1)
    CFITSIO.fits_read_col(col.hdu.fitsfile, col.properties.colnum, i, 1, A)
    return only(A)
end

function Base.getindex(col::TableHDUColumn{<:AbstractArray}, i::Int)
    ALLOW_SCALAR[] || error("Scalar access to FITS columns is disabled. Enable it by calling `allow_scalar()`.")
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
