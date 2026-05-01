struct GroupedHDU <: FITSIO.HDU
    fitsfile::FITSIO.FITSFile
    ext::Int
end

function Base.read(hdu::GroupedHDU)
    FITSIO.assert_open(hdu)
    FITSIO.fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = FITSIO.fits_get_img_size(hdu.fitsfile) |> Tuple
    @assert first(sz) == 0
    h = FITSIO.read_header(hdu)
    ngroups = h["GCOUNT"]
    pcount = h["PCOUNT"]

    pnames = map(1:pcount) do i
        h["PTYPE$i"]
    end
    for i in 1:length(pnames)
        while pnames[i] in pnames[1:i-1]
            pnames[i] = "_" * pnames[i]
        end
    end

    ps_scalezero = ntuple(pcount) do i
        get(h, "PSCAL$i", 1.0), get(h, "PZERO$i", 0.0)
    end
    haskey(h, "BSZERO") && @warn "FITS header contains BSZERO keyword, probably a typo"

    ptype = NTuple{pcount, Float32}
    result_grp = NamedTuple{Tuple(Symbol.(pnames))}(ntuple(_ -> Float64[], pcount))
    result_data = Array{Float32}(undef, Base.tail(sz)..., ngroups)

    T = FITSIO.type_from_bitpix(FITSIO.fits_get_img_equivtype(hdu.fitsfile))
    @assert T == Float32

    buf_grp = Array{Cfloat}(undef, pcount)
    buf_data = Array{Cfloat}(undef, Base.tail(sz))

    _fill_data!(ngroups, hdu, buf_grp, buf_data, result_grp, result_data, ptype, ps_scalezero)

    return (; result_grp..., DATA=result_data)
end

function _fill_data!(ngroups, hdu, buf_grp, buf_data, result_grp, result_data, ::Type{ptype}, ps_scalezero) where {ptype}
    for groupix in 1:ngroups
        status = Ref{Cint}(0)
        @ccall FITSIO.libcfitsio.ffggpe(
            hdu.fitsfile.ptr::Ptr{Cvoid},
            groupix::Clong,
            1::Clong,  # firstelem
            length(buf_grp)::Clong,
            buf_grp::Ptr{Cfloat},
            status::Ref{Cint},
        )::Cint
        FITSIO.fits_assert_ok(status[])

        anynul = Ref{Cint}(0)
        @ccall FITSIO.libcfitsio.ffgpve(
            hdu.fitsfile.ptr::Ptr{Cvoid},
            groupix::Clong,
            1::Clong,  # firstelem
            length(buf_data)::Clong,
            0.0::Cfloat,  # nulval
            buf_data::Ptr{Cfloat},
            anynul::Ref{Cint},
            status::Ref{Cint},
        )::Cint
        @assert anynul[] == 0
        FITSIO.fits_assert_ok(status[])

        map(Tuple(result_grp), ptype(buf_grp), ps_scalezero) do vres, val, scalezero
            push!(vres, muladd(val, scalezero...))
        end
        selectdim(result_data, ndims(result_data), groupix) .= buf_data
    end
end

# --- Lazy random-groups reader (mirrors fast_column_read.jl for BINTABLEs) ---

using Mmap

"""
    MmapGroupContext

Mmap'd handle to a random-groups data unit. `data` is the entire file's bytes;
`data_start` is the byte offset (1-based, suitable for `pointer(data, i)`) of
the first group's first parameter; `group_bytes` is the per-group stride in
bytes (params + data).
"""
struct MmapGroupContext
    filepath::String
    data::Vector{UInt8}
    data_start::Int
    group_bytes::Int
    pcount::Int
    data_elem_size::Int       # |BITPIX|/8 (typically 4 for BITPIX=-32)
    data_nelem_per_group::Int # prod(NAXIS2..NAXISm)
    ngroups::Int
end

function _mmap_groupedhdu_context(hdu::GroupedHDU)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    filepath = fits_file_name(hdu.fitsfile)
    isfile(filepath) || return nothing

    sz = FITSIO.fits_get_img_size(hdu.fitsfile) |> Tuple
    @assert first(sz) == 0 "Not a random-groups HDU (NAXIS1 != 0)"
    h = FITSIO.read_header(hdu)
    bitpix = h["BITPIX"]
    elem_size = abs(bitpix) ÷ 8
    @assert bitpix == -32 "Only BITPIX=-32 (Float32) random groups supported; got BITPIX=$bitpix"

    pcount = h["PCOUNT"]
    ngroups = h["GCOUNT"]
    data_nelem = prod(Base.tail(sz))   # NAXIS2 × NAXIS3 × ... × NAXISm
    group_bytes = (pcount + data_nelem) * elem_size

    addr = _fits_get_hduaddr(hdu.fitsfile)
    io = open(filepath, "r")
    data = Mmap.mmap(io, Vector{UInt8})
    close(io)
    ccall(:madvise, Cint, (Ptr{UInt8}, Csize_t, Cint), pointer(data), length(data), 2)  # MADV_SEQUENTIAL

    return MmapGroupContext(filepath, data, addr.datastart, group_bytes,
                            pcount, elem_size, Int(data_nelem), Int(ngroups))
end

"""
    MmapGroupDataRow{T}

A lazy `AbstractVector{T}` view of one random-group's data block, flat (the
NAXIS2..NAXISm dims are not exposed at this layer — consumers reshape).
Carries its own `pos` (byte offset of this group's data start) and `repeat`
(element count). Mirrors `MmapArrayRow` from fast_column_read.jl.
"""
struct MmapGroupDataRow{T<:AbstractFloat} <: AbstractVector{T}
    ctx::MmapGroupContext
    pos::Int
    repeat::Int
end

Base.size(r::MmapGroupDataRow) = (r.repeat,)
Base.IndexStyle(::Type{<:MmapGroupDataRow}) = IndexLinear()

@inline function Base.getindex(r::MmapGroupDataRow{T}, i::Int) where {T<:AbstractFloat}
    @boundscheck checkbounds(r, i)
    ctx = r.ctx
    off = r.pos + (i - 1) * ctx.data_elem_size
    raw = GC.@preserve ctx ntoh(unsafe_load(Ptr{Float32}(pointer(ctx.data, off))))
    return T(raw)
end

"""
    MmapGroupColumn{ET}

A lazy `AbstractVector{ET}` over all groups of a random-groups HDU. Unifies:

- **scalar PTYPE columns** (`ET <: AbstractFloat`): one Float32 per group at
  `byte_offset_in_group`, scaled by `physical = pzero + pscal * raw` (FITS
  Standard §6.1.2). The output type `ET` is picked at construction time by
  `_column_eltype(pscal, pzero)` — `Float32` if PSCAL/PZERO are trivially
  1.0/0.0 (the on-disk Float32 *is* the physical value, no precision lost),
  otherwise `Float64` (muladd in Float32 would round PZERO to 24 bits —
  fatal for e.g. DATE's JD-at-midnight).

- **the DATA array column** (`ET = MmapGroupDataRow{Float32}`): `getindex(g)`
  returns a `MmapGroupDataRow{Float32}` view of group `g`'s data block,
  flat-indexed. No scaling.

`byte_offset_in_group` is the column's offset from the start of each group:
`(i-1) * data_elem_size` for param `i`, `pcount * data_elem_size` for DATA.
`repeat` is 1 for scalar params, `data_nelem_per_group` for DATA.

Mirrors `MmapColumn{ET, T}` from fast_column_read.jl — same unification
trick, two `getindex` methods on `ET`.
"""
struct MmapGroupColumn{ET} <: AbstractVector{ET}
    ctx::MmapGroupContext
    byte_offset_in_group::Int
    repeat::Int
    pscal::Float64
    pzero::Float64
end

Base.size(c::MmapGroupColumn) = (c.ctx.ngroups,)
Base.IndexStyle(::Type{<:MmapGroupColumn}) = IndexLinear()

# Pick the smallest lossless output type for a scalar PTYPE column.
_column_eltype(pscal::Float64, pzero::Float64) =
    (pscal == 1.0 && pzero == 0.0) ? Float32 : Float64

# Scalar, no scaling: bit-exact return of the on-disk Float32.
@inline function Base.getindex(c::MmapGroupColumn{Float32}, g::Int)
    @boundscheck checkbounds(c, g)
    ctx = c.ctx
    pos = ctx.data_start + (g - 1) * ctx.group_bytes + c.byte_offset_in_group + 1
    GC.@preserve ctx ntoh(unsafe_load(Ptr{Float32}(pointer(ctx.data, pos))))
end

# Scalar, with scaling: muladd promotes Float32 raw to Float64 automatically;
# returning Float64 preserves precision when PZERO has more than ~24 bits of
# significand (canonical case: DATE).
@inline function Base.getindex(c::MmapGroupColumn{Float64}, g::Int)
    @boundscheck checkbounds(c, g)
    ctx = c.ctx
    pos = ctx.data_start + (g - 1) * ctx.group_bytes + c.byte_offset_in_group + 1
    raw = GC.@preserve ctx ntoh(unsafe_load(Ptr{Float32}(pointer(ctx.data, pos))))
    return muladd(raw, c.pscal, c.pzero)
end

# Array: return a per-group row view.
@inline function Base.getindex(c::MmapGroupColumn{MmapGroupDataRow{T}}, g::Int) where {T<:AbstractFloat}
    @boundscheck checkbounds(c, g)
    ctx = c.ctx
    pos = ctx.data_start + (g - 1) * ctx.group_bytes + c.byte_offset_in_group + 1
    return MmapGroupDataRow{T}(ctx, pos, c.repeat)
end

"""
    lazycolumntable(hdu::GroupedHDU) → NamedTuple

Lazy/mmap'd analogue of `Base.read(hdu)` and of `lazycolumntable(::TableHDU)`.
Returns a NamedTuple keyed by PTYPE names (with `_`-prefix deduplication
matching the eager reader), holding `MmapGroupColumn{ET}` columns. Each
scalar column's `ET` is auto-selected: `Float32` if PSCAL/PZERO are 1.0/0.0,
else `Float64` (see `_column_eltype`). The `:DATA` column has
`ET = MmapGroupDataRow{Float32}`.

Falls back to `nothing` (caller's responsibility to detect) if the file can't
be mmapped.
"""
function lazycolumntable(hdu::GroupedHDU)
    ctx = _mmap_groupedhdu_context(hdu)
    ctx === nothing && return nothing
    h = FITSIO.read_header(hdu)

    pcount = ctx.pcount
    pnames = map(1:pcount) do i
        h["PTYPE$i"]
    end
    # Same dedup as eager Base.read (lines 18-22 of this file).
    for i in 1:pcount
        while pnames[i] in pnames[1:i-1]
            pnames[i] = "_" * pnames[i]
        end
    end

    params = ntuple(pcount) do i
        pscal = Float64(get(h, "PSCAL$i", 1.0))
        pzero = Float64(get(h, "PZERO$i", 0.0))
        offset = (i - 1) * ctx.data_elem_size
        ET = _column_eltype(pscal, pzero)
        MmapGroupColumn{ET}(ctx, offset, 1, pscal, pzero)
    end

    data_col = MmapGroupColumn{MmapGroupDataRow{Float32}}(
        ctx, ctx.pcount * ctx.data_elem_size, ctx.data_nelem_per_group, 1.0, 0.0)

    sym_pnames = Tuple(Symbol.(pnames))
    return NamedTuple{(sym_pnames..., :DATA)}((params..., data_col))
end
