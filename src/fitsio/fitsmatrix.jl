using StaticArrays, AxisKeys
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
