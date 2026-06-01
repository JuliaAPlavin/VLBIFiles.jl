using Mmap
using FITSIO.Libcfitsio: fits_get_coltype, fits_get_colnum, fits_get_num_rows, fits_file_name, fits_movabs_hdu

struct MmapTableContext
    filepath::String
    data::Vector{UInt8}
    data_start::Int
    row_bytes::Int
    nrows::Int
end

struct MmapColumn{ET, T<:Number} <: AbstractVector{ET}
    ctx::MmapTableContext
    col_offset::Int
    repeat::Int
end

struct MmapArrayRow{T<:Number} <: AbstractVector{T}
    ctx::MmapTableContext
    pos::Int
    repeat::Int
end

function _mmap_column(ctx::MmapTableContext, hdu::TableHDU, colnum::Integer)
    typecode, repeat, width = fits_get_coltype(hdu.fitsfile, colnum)
    col_offset = _compute_column_byte_offset(hdu.fitsfile, colnum)
    T = get(_CFITSIO_TYPECODE_TO_JULIA, typecode, nothing)
    T === nothing && return nothing  # unsupported type (e.g. string)
    if repeat == 1
        return MmapColumn{T, T}(ctx, col_offset, 1)
    else
        return MmapColumn{MmapArrayRow{T}, T}(ctx, col_offset, repeat)
    end
end

Base.size(col::MmapColumn) = (col.ctx.nrows,)
Base.size(row::MmapArrayRow) = (row.repeat,)
Base.IndexStyle(::Type{<:MmapArrayRow}) = IndexLinear()

# --- Scalar getindex ---

@inline function Base.getindex(col::MmapColumn{T, T}, i::Int) where {T<:Number}
    @boundscheck checkbounds(col, i)
    pos = col.ctx.data_start + col.col_offset + (i - 1) * col.ctx.row_bytes + 1
    GC.@preserve col ntoh(unsafe_load(Ptr{T}(pointer(col.ctx.data, pos))))
end

# --- Array getindex ---

@inline function Base.getindex(row::MmapArrayRow{T}, i::Int) where {T<:Number}
    @boundscheck checkbounds(row, i)
    GC.@preserve row ntoh(unsafe_load(Ptr{T}(pointer(row.ctx.data, row.pos + (i - 1) * sizeof(T)))))
end

@inline function Base.getindex(col::MmapColumn{MmapArrayRow{T}, T}, i::Int) where {T<:Number}
    @boundscheck checkbounds(col, i)
    pos = col.ctx.data_start + col.col_offset + (i - 1) * col.ctx.row_bytes + 1
    MmapArrayRow{T}(col.ctx, pos, col.repeat)
end

# --- Fast bulk collect via single-pass pread ---

Base.collect(col::MmapColumn) = only(materialize_columns((col,)))
Base.copyto!(dest::AbstractVector, col::MmapColumn) = (_multi_collect((col,), col.ctx, (dest,)); dest)

# --- Multi-column single-pass collect ---

function materialize_columns(cols::Union{Tuple, NamedTuple})
    @assert all(c -> c isa MmapColumn, cols)
    @assert allequal(c.ctx for c in cols)
    _multi_collect(cols, first(cols).ctx)
end

_materialized_similar(col::MmapColumn{T, T}) where {T<:Number} = Vector{T}(undef, length(col))
_materialized_similar(col::MmapColumn{MmapArrayRow{T}, T}) where {T<:Number} = Vector{Vector{T}}(undef, length(col))

function _multi_collect(cols, ctx::MmapTableContext, dests=map(_materialized_similar, cols))
    nrows = ctx.nrows
    row_bytes = ctx.row_bytes
    bufsize = 32 * 1024 * 1024
    rows_per_buf = max(1, bufsize ÷ row_bytes)
    buf = Vector{UInt8}(undef, rows_per_buf * row_bytes)

    fd = ccall(:open, Cint, (Cstring, Cint), ctx.filepath, 0)
    fd == -1 && error("Failed to open $(ctx.filepath)")
    try
        dest_idx = 1
        rows_remaining = nrows
        file_pos = Int64(ctx.data_start)
        while rows_remaining > 0
            rows_this = min(rows_per_buf, rows_remaining)
            nbytes = ccall(:pread, Cssize_t, (Cint, Ptr{UInt8}, Csize_t, Int64),
                fd, buf, rows_this * row_bytes, file_pos)
            nbytes < rows_this * row_bytes && error("Short read")
            bp = pointer(buf)
            for r in 0:(rows_this - 1)
                row_ptr = bp + r * row_bytes
                _extract_row!(dests, cols, row_ptr, dest_idx)
                dest_idx += 1
            end
            file_pos += rows_this * row_bytes
            rows_remaining -= rows_this
        end
    finally
        ccall(:close, Cint, (Cint,), fd)
    end
    return dests
end

@inline function _extract_row!(dests, cols, row_ptr::Ptr{UInt8}, idx::Int)
    @inline map(dests, cols) do dest, col
        _extract_one!(dest, col, row_ptr, idx)
    end
end

@inline function _extract_one!(dest::AbstractVector, col::MmapColumn{T, T}, row_ptr::Ptr{UInt8}, idx::Int) where {T<:Number}
    @inbounds dest[idx] = ntoh(unsafe_load(Ptr{T}(row_ptr + col.col_offset)))
end

@inline function _extract_one!(dest::AbstractVector{Vector{T}}, col::MmapColumn{MmapArrayRow{T}, T}, row_ptr::Ptr{UInt8}, idx::Int) where {T<:Number}
    v = Vector{T}(undef, col.repeat)
    pp = row_ptr + col.col_offset
    @inbounds for j in 1:col.repeat
        v[j] = ntoh(unsafe_load(Ptr{T}(pp + (j - 1) * sizeof(T))))
    end
    @inbounds dest[idx] = v
end

@inline function _extract_one!(dest::AbstractVector, col::MmapColumn{MmapArrayRow{T}, T}, row_ptr::Ptr{UInt8}, idx::Int) where {T<:Number}
    @inbounds dest[idx] = col[idx]
end

# --- Helpers ---

function _fits_get_hduaddr(fitsfile)
    headstart = Ref{Clonglong}(0)
    datastart = Ref{Clonglong}(0)
    dataend = Ref{Clonglong}(0)
    status = Ref{Cint}(0)
    ccall((:ffghadll, FITSIO.libcfitsio), Cint,
        (Ptr{Cvoid}, Ref{Clonglong}, Ref{Clonglong}, Ref{Clonglong}, Ref{Cint}),
        fitsfile.ptr, headstart, datastart, dataend, status)
    status[] != 0 && error("ffghadll failed with status $(status[])")
    return (headstart=Int(headstart[]), datastart=Int(datastart[]), dataend=Int(dataend[]))
end

function _column_byte_width(fitsfile, colnum::Integer)
    typecode, repeat, width = fits_get_coltype(fitsfile, colnum)
    # For string columns (typecode=16), width equals repeat (display width),
    # but actual storage is 1 byte per character = repeat bytes total.
    typecode == 16 ? repeat : repeat * width
end

_compute_column_byte_offset(fitsfile, target_colnum::Integer) =
    sum(c -> _column_byte_width(fitsfile, c), 1:(target_colnum - 1); init=0)

const _CFITSIO_TYPECODE_TO_JULIA = Dict{Int, DataType}(
    # TLOGICAL (14) deliberately excluded: FITS stores 'T'/'F'/'\0' bytes,
    # which cannot be correctly interpreted via ntoh(unsafe_load(Ptr{Bool}(...)))
    11 => UInt8,    # TBYTE
    12 => Int16,    # TSHORT
    21 => Int32,    # TLONG (equivalent)
    41 => Int32,    # TLONG
    81 => Int64,    # TLONGLONG
    42 => Float32,  # TFLOAT
    82 => Float64,  # TDOUBLE
)

function _mmap_table_context(hdu::TableHDU)
    _ensure_hdu_active!(hdu)
    filepath = fits_file_name(hdu.fitsfile)
    isfile(filepath) || return nothing

    addr = _fits_get_hduaddr(hdu.fitsfile)
    fhead = FITSIO.read_header(hdu)
    row_bytes = fhead["NAXIS1"]
    nrows = fits_get_num_rows(hdu.fitsfile)

    # Validate that our column byte width computation agrees with NAXIS1
    ncols = fhead["TFIELDS"]
    computed_row_bytes = sum(_column_byte_width(hdu.fitsfile, c) for c in 1:ncols)
    computed_row_bytes == row_bytes || error("computed column byte widths disagree with NAXIS1")

    io = open(filepath, "r")
    data = Mmap.mmap(io, Vector{UInt8})
    close(io)

    # ccall(:madvise, Cint, (Ptr{UInt8}, Csize_t, Cint), pointer(data), length(data), 1)  # MADV_RANDOM=1
    ccall(:madvise, Cint, (Ptr{UInt8}, Csize_t, Cint), pointer(data), length(data), 2) # MADV_SEQUENTIAL=2

    return MmapTableContext(filepath, data, addr.datastart, row_bytes, nrows)
end
