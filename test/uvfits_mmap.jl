@testitem "UVFITS lazy random-groups reader vs eager" begin
    using VLBIFiles: FITSIO

    fixtures = [
        "0332-391.uvfits",
        "SR1_3C279_2017_101_hi_hops_netcal_StokesI.uvfits",
        "hops_3600_OJ287_LO+HI.medcal_dcal_full.uvfits",
        "datafile_01-01_230GHz.uvfits",
        "mwa_1061316296.uvfits",
        "paper_zen.uvfits",
    ]

    for fname in fixtures
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

            # Each parameter column: same length, same values element-wise. Eager is always
            # Vector{Float64}; lazy is the honest on-disk precision — Float32 for a pure-rescale
            # column (UU/VV/WW: PZERO=0), Float64 only when PZERO needs the headroom (DATE). So
            # compare eager rounded to the lazy column's precision.
            for k in keys(eager)
                k === :DATA && continue
                @test isequal(eltype(lazy[k]).(eager[k]), lazy[k])
            end

            # DATA: eager is `Array{Float32, 7}` with groups along the last axis;
            # lazy is a column whose elements are per-group N-D KeyedArrays
            # (COMPLEX, STOKES, FREQ, IF, RA, DEC). Values are column-major identical,
            # so compare per-group modulo shape via `vec` — at head, tail, and middle.
            edata = eager[:DATA]
            ldata_col = lazy[:DATA]
            @test length(ldata_col) == size(edata, ndims(edata))
            for g in [1, length(ldata_col), length(ldata_col) ÷ 2 + 1]
                e_g = vec(selectdim(edata, ndims(edata), g))
                l_g = ldata_col[g]
                @test isequal(e_g, vec(l_g))
            end
        end
    end
end
