using Test
using VLBIFiles
using VLBIFiles: FITSIO

UVFITS_FIXTURES = [
    "0332-391.uvfits",
    "SR1_3C279_2017_101_hi_hops_netcal_StokesI.uvfits",
    "hops_3600_OJ287_LO+HI.medcal_dcal_full.uvfits",
    "datafile_01-01_230GHz.uvfits",
    "mwa_1061316296.uvfits",
    "paper_zen.uvfits",
]

@testset "UVFITS lazy random-groups reader vs eager" begin
    for fname in UVFITS_FIXTURES
        path = joinpath(pkgdir(VLBIFiles), "test", "data", fname)
        @testset "$fname" begin
            # Eager path (existing)
            eager = FITSIO.FITS(path) do f
                read(VLBIFiles.GroupedHDU(f.fitsfile, 1))
            end

            # Lazy path (new)
            lazy = FITSIO.FITS(path) do f
                VLBIFiles.lazycolumntable(VLBIFiles.GroupedHDU(f.fitsfile, 1))
            end

            # Same parameter columns in the same order — both readers build
            # the NamedTuple in PTYPE order (PTYPE1..PCOUNT) followed by :DATA,
            # so we compare ordered. Set-based equality would silently mask
            # reorder bugs.
            @test keys(eager) == keys(lazy)

            # Each parameter column: same length, same values element-wise.
            # `isequal` on two AbstractVectors compares pairwise without
            # requiring the same eltype — eager is Vector{Float64}, lazy may
            # be MmapGroupColumn{Float32} or {Float64} depending on PSCAL/PZERO.
            # For the no-scaling case, Float64(Float32_raw) == Float32_raw is
            # bit-exact, so the cross-type comparison still succeeds.
            for k in keys(eager)
                k === :DATA && continue
                @test isequal(eager[k], lazy[k])
            end

            # DATA: eager is `Array{Float32, 7}` with groups along the last
            # axis; lazy is `MmapGroupColumn{MmapGroupDataRow{Float32}}` of
            # length GCOUNT whose elements are flat per-group views. Compare
            # per-group at the head, tail, and middle (no collect needed —
            # `vec(selectdim(...))` is a view, `isequal` walks both lazily).
            edata = eager[:DATA]
            ldata_col = lazy[:DATA]
            @test length(ldata_col) == size(edata, ndims(edata))
            for g in [1, length(ldata_col), length(ldata_col) ÷ 2 + 1]
                e_g = vec(selectdim(edata, ndims(edata), g))
                l_g = ldata_col[g]
                @test isequal(e_g, l_g)
            end
        end
    end
end
