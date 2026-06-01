using NamedDims: NamedDimsArray

# One row's visibilities as a lazy (stokes, channel) matrix. The flat data block is reshaped on
# construction into a *named* (COMPLEX, stokes, channel) array — COMPLEX (re, im[, wt]) varies fastest —
# so re/im read by name: block[COMPLEX=1] / block[COMPLEX=2]. reshape of an mmap-backed vector is a
# zero-copy ReshapedArray; NamedDimsArray adds only compile-time dim names (no runtime cost, verified
# 0-alloc and identical speed to positional indexing); both throw if the length isn't n_complex·n_stokes·n_chan.
const _BLOCK_AXES = (:COMPLEX, :stokes, :channel)
struct FitsVMatrix{T,P} <: AbstractMatrix{Complex{T}}
    block::P
end
FitsVMatrix(flat::AbstractArray, n_complex, n_stokes, n_chan) =
    FitsVMatrix(NamedDimsArray(reshape(flat, n_complex, n_stokes, n_chan), _BLOCK_AXES))
FitsVMatrix(block::AbstractArray{T,3}) where T = FitsVMatrix{T,typeof(block)}(block)
Base.size(v::FitsVMatrix) = (size(v.block, :stokes), size(v.block, :channel))
Base.@propagate_inbounds Base.getindex(v::FitsVMatrix, s::Int, k::Int) =
    complex(v.block[COMPLEX=1, stokes=s, channel=k], v.block[COMPLEX=2, stokes=s, channel=k])

# Embedded weight (COMPLEX axis size 3, used by UVFITS and FITS-IDI alike): the 3rd COMPLEX entry of
# the named (COMPLEX, stokes, channel) block — one weight per (stokes, channel). Same _BLOCK_AXES as
# FitsVMatrix, read by name.
struct FitsWMatrixEmb{T,P} <: AbstractMatrix{T}
    block::P
end
FitsWMatrixEmb(flat::AbstractArray, n_complex, n_stokes, n_chan) =
    FitsWMatrixEmb(NamedDimsArray(reshape(flat, n_complex, n_stokes, n_chan), _BLOCK_AXES))
FitsWMatrixEmb(block::AbstractArray{T,3}) where T = FitsWMatrixEmb{T,typeof(block)}(block)
Base.size(w::FitsWMatrixEmb) = (size(w.block, :stokes), size(w.block, :channel))
Base.@propagate_inbounds Base.getindex(w::FitsWMatrixEmb, s::Int, k::Int) = w.block[COMPLEX=3, stokes=s, channel=k]

# FITS-IDI WEIGHT column: one value per (stokes, band), reshaped to a named (stokes, band) array and
# broadcast across the band's channels. reshape throws if the block length isn't n_stokes·n_IF.
const _WSEP_AXES = (:stokes, :band)
struct FitsWMatrixSep{T,P} <: AbstractMatrix{T}
    block::P
    nchan_per_IF::Int
    n_chan::Int
end
FitsWMatrixSep(flat::AbstractArray, n_stokes, nchan_per_IF, n_chan) =
    FitsWMatrixSep(NamedDimsArray(reshape(flat, n_stokes, n_chan ÷ nchan_per_IF), _WSEP_AXES), nchan_per_IF, n_chan)
FitsWMatrixSep(block::AbstractArray{T,2}, nchan_per_IF, n_chan) where T = FitsWMatrixSep{T,typeof(block)}(block, nchan_per_IF, n_chan)
Base.size(w::FitsWMatrixSep) = (size(w.block, :stokes), w.n_chan)
Base.@propagate_inbounds Base.getindex(w::FitsWMatrixSep, s::Int, k::Int) =
    w.block[stokes=s, band=(k-1) ÷ w.nchan_per_IF + 1]


# Lazy per-channel 2×2 over a (stokes,channel) matrix, using the feed→position layout from the matrix's
# stokes axis. `absent` fills positions whose stokes is missing; `value` post-processes present
# entries (identity for visibilities; the per-feed weight policy for weights). Both are folded into
# getindex, so each access is one allocation-free SMatrix.
struct CoherencyView{T,C,F} <: AbstractVector{SMatrix{2,2,T,4}}
    mat::C
    slot::SMatrix{2,2,Int8,4}
    n_chan::Int
    absent::T
    value::F
end
Base.size(v::CoherencyView) = (v.n_chan,)
Base.@propagate_inbounds Base.getindex(v::CoherencyView, k::Int) =
    map(i -> i == 0 ? v.absent : v.value(v.mat[i, k]), v.slot)
# Specialized force-inlined iterate (`@propagate_inbounds` ⇒ `@inline`): the generic AbstractArray
# iterate is cost-gated by the deep lazy getindex chain and won't inline into a loop, leaving the
# loop-invariant pointer/stride setup un-hoisted (2-3× slower). This makes `zip` over the lazy
# visibility/weight columns as fast as indexed access. See Fringy gridding (`fringe_grids`).
Base.@propagate_inbounds Base.iterate(v::CoherencyView, c::Int=1) =
    c > v.n_chan ? nothing : (@inbounds(v[c]), c + 1)
# A missing complex visibility is NaN+NaNim (both parts, so isnan holds and nothing inspects a stray 0im).
_absent(::Type{Complex{T}}) where T = Complex{T}(NaN, NaN)
_absent(::Type{T}) where T = T(NaN)
# Two entry points: a slot-taking core (a pure view over a precomputed slot), and a column method that
# computes the slot once from the shared stokes axis (identical across a file's records) and reuses it.
@inline coherencymatrices(mat, slot::SMatrix{2,2,Int8}; absent=_absent(eltype(mat)), value=identity) =
    CoherencyView(mat, slot, size(mat, 2), absent, value)
function coherencymatrices(col::AbstractVector; absent=_absent(eltype(eltype(col))), value=identity)
    slot = VLBIData.feed_slot_map(axiskeys(first(col), :stokes))   # stokes axis is shared across the column
    mapview(c -> coherencymatrices(c, slot; absent, value), col)
end
