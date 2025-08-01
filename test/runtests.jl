using TestItems
using TestItemRunner
@run_package_tests


@testitem "generic loading" begin
    cd(dirname(@__FILE__))
    @test VLBI.guess_type("./data/map.fits") == VLBI.FitsImage
    @test VLBI.guess_type("./data/vis.fits") == VLBI.UVData
    @test VLBI.guess_type("./data/vis_multichan.vis") == VLBI.UVData
    @test VLBI.guess_type("./data/difmap_model.mod") == MultiComponentModel
    @test VLBI.guess_type("./data/difmap_model_empty.mod") == MultiComponentModel
end

@testitem "img don't read data" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    cd(dirname(@__FILE__))

    img = VLBI.load("./data/map.fits", read_data=false)
    @test_throws "not loaded" KeyedArray(img)

    @test VLBI.pixel_size(img) ≈ 0.2u"mas"  rtol=1e-5
    @test VLBI.pixel_steps(img) ≈ [-0.2u"mas", 0.2u"mas"]  rtol=1e-5
    @test VLBI.pixel_area(img) ≈ 0.04u"mas^2"  rtol=1e-5
    @test frequency(img) ≈ 4.344u"GHz"

    bm = beam(img)
    @test VLBI.load(Beam, "./data/map.fits") === bm
    @test effective_area(bm) ≈ 6.5519451u"mas^2"  rtol=1e-5
    @test fwhm_max(bm) ≈ 4.02887356u"mas"  rtol=1e-5
    @test fwhm_min(bm) ≈ 1.43523228u"mas"  rtol=1e-5
    @test position_angle(bm) ≈ -0.046499863  rtol=1e-5
    @test intensity(ModelComponent(bm))(SVector(0.5u"mas", 0u"mas")) ≈ 0.714721  rtol=1e-5
end
    
@testitem "img read data" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    cd(dirname(@__FILE__))

    img = VLBI.load("./data/map.fits", read_data=true)
    @test img.data == VLBI.load(VLBI.FitsImage, "./data/map.fits").data
    KA = KeyedArray(img)
    @test KA == VLBI.load(KeyedArray, "./data/map.fits")
    @test size(KA) == (512, 512)
    @test mean(KA) ≈ 2.4069064f-5
    @test maximum(KA) ≈ 0.021357307f0
    @test KA[123, 456] ≈ -7.372097f-5
    @test axiskeys(KA, :ra)  .|> ustrip ≈ 51:-0.2:-51.2  atol=1e-3
    @test axiskeys(KA, :dec) .|> ustrip ≈ -51.2:0.2:51   atol=1e-3
    @test axiskeys(KA, :ra) isa AbstractRange
    @test axiskeys(KA, :dec) isa AbstractRange
end

@testitem "img save" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using AxisKeys

    tmpf = tempname()
    
    img = KeyedArray(rand(512, 256), ra=range(-10, 10, length=512)u"mas", dec=range(-20, 20, length=256)u"mas")
    @test_throws Exception VLBI.save(tmpf, img)

    img = KeyedArray(rand(512, 256), ra=range(-10, 10, length=512)u"mas", dec=range(-20, 20, length=256)u"mas")*u"Jy"
    VLBI.save(tmpf, img; freq=123)

    img_r = VLBI.load(tmpf).data
    @test axiskeys(img_r, 1) ≈ axiskeys(img, 1)
    @test axiskeys(img_r, 2) ≈ axiskeys(img, 2)
    @test AxisKeys.keyless_unname(img_r) ≈ AxisKeys.keyless_unname(ustrip.(u"Jy", img))
    @test_broken VLBI.frequency(VLBI.load(tmpf)) == 123

    img = KeyedArray(rand(512, 256), ra=range(-5, 10, length=512)u"mas", dec=range(-7, 1, length=256)u"mas")*u"Jy"
    VLBI.save(tmpf, img)

    img_r = VLBI.load(tmpf).data
    @test axiskeys(img_r, 1) ≈ axiskeys(img, 1)
    @test axiskeys(img_r, 2) ≈ axiskeys(img, 2)
    @test AxisKeys.keyless_unname(img_r) ≈ AxisKeys.keyless_unname(ustrip.(u"Jy", img))
end

@testitem "img read clean" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    cd(dirname(@__FILE__))

    mod = VLBI.load(MultiComponentModel, "./data/map.fits")
    @test length(components(mod)) == 361
    @test sum(flux, components(mod)) ≈ 0.038659f0u"Jy"
    @test mean(first ∘ coords, components(mod)) ≈ -0.913573508654587u"mas"
    @test mean(last ∘ coords, components(mod)) ≈ 8.145706523588233u"mas"
end

@testitem "img stacked" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    cd(dirname(@__FILE__))

    img = VLBI.load("./data/map_stacked.fits", read_data=true)
    KA = KeyedArray(img)
    @test size(KA) == (512, 512)
    @test maximum(KA) ≈ 0.5073718945185344
    @test KA[123, 456] ≈ -1.0848011697817128e-5
    @test axiskeys(KA, :ra)  .|> ustrip ≈ 25.5:-0.1:-25.6  atol=1e-3
    bm = beam(img)
    @test fwhm_min(bm) == fwhm_max(bm)

    @test_broken mod = VLBI.load(MultiComponentModel, "./data/map_stacked.fits")
end

@testitem "img clean/residual/combined" begin
    using RectiGrids, AxisKeysExtra, Interpolations
    using Unitful, UnitfulAngles, UnitfulAstro
    using StaticArrays

    img = VLBI.load("./data/map.fits")
    nativebeam = Beam(img)
    newbeam = Beam(CircularGaussian, σ = 5u"mas")
    @test with_axiskeys(grid)(img) |> size == (512, 512)
    
    @test VLBI.image_stored(KeyedArray, img) === img.data
    @test VLBI.image_stored(KeyedArray, img)(0, 0) ≈ 0.021357307f0

    @test VLBI.image_clean(img; beam=nativebeam)(SVector(0,0)) ≈ 0.021302195437816945
    @test VLBI.image_clean(img; beam=nativebeam)(SVector(3,0)u"mas") ≈ 0.0003380008663429093
    @test VLBI.image_clean(img; beam=newbeam)(SVector(3,0)u"mas") ≈ 0.02823164807298988
    xy = rand(SVector{2,Float64})u"mas"
    @test VLBI.image_clean(img)(xy) == VLBI.image_clean(img; beam=nativebeam)(xy)

    @test VLBI.image_clean_with_residual(img)(xy) == VLBI.image_clean_with_residual(img; beam=nativebeam)(xy)
    @test VLBI.image_clean_with_residual(KeyedArray, img) == VLBI.image_clean(KeyedArray, img) + VLBI.image_residual(KeyedArray, img)
    @test VLBI.image_clean_with_residual(KeyedArray, img; beam=newbeam) == VLBI.image_clean(KeyedArray, img; beam=newbeam) + VLBI.image_residual(KeyedArray, img)

    xy = with_axiskeys(grid)(SVector, img) |> rand
    @testset "lazy-mat" for (lazy, materialized) in Any[
        (VLBI.image_stored(img), VLBI.image_stored(KeyedArray, img)),
        (VLBI.image_clean(img), VLBI.image_clean(KeyedArray, img)),
        (VLBI.image_clean(img; beam=newbeam), VLBI.image_clean(KeyedArray, img; beam=newbeam)),
        (VLBI.image_residual(img), VLBI.image_residual(KeyedArray, img)),
        (VLBI.image_clean_with_residual(img), VLBI.image_clean_with_residual(KeyedArray, img)),
        (VLBI.image_clean_with_residual(img; beam=newbeam), VLBI.image_clean_with_residual(KeyedArray, img; beam=newbeam)),
    ]
        @test lazy(SVector(0, 0)) > 0
        @test lazy(SVector(0, 0)) ≈ materialized(0, 0)
        @test lazy(xy) ≈ materialized(xy...)
    end
end

@testitem "img nonstandard header names" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    cd(dirname(@__FILE__))

    img = VLBI.load("./data/sampling_mean.fits", read_data=true)
    KA = KeyedArray(img)
    @test size(KA) == (256, 256)
    @test maximum(KA) ≈ 0.0022970616631588308
    @test KA[123, 222] ≈ 1.0227093348779382e-6
    @test axiskeys(KA, :ra)  .|> ustrip ≈ -0.3984375:0.003125:0.3984375  atol=1e-3
    @test_throws "key \"BMAJ\" not found" beam(img)
end

@testitem "uvf simple" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIFiles.Uncertain
    using Statistics
    using StaticArrays
    using DataManipulation
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/vis.fits")
    @test uv.header.object == "J1033+6051"
    @test Date(uv.header) === Date(2010, 12, 24)
    @test uv.header.stokes == [:RR, :LL, :RL, :LR]
    @test frequency(uv.header) ≈ 15.33u"GHz" rtol=1e-2
    @test length(uv.freq_windows) == 8
    @test length(uv.ant_arrays) == 1
    antarr = only(uv.ant_arrays)
    @test antarr.name == "VLBA"
    @test map(a -> a.name, antarr.antennas) |> collect == [:BR, :FD, :HN, :KP, :LA, :MK, :NL, :OV, :PT, :SC]

    @test frequency(uv.freq_windows[1]) == 1.5329522f10u"Hz"
    @test frequency(uv.freq_windows[1], VLBIFiles.Interval) == VLBIFiles.Interval(1.5329522f10u"Hz", 1.5337521f10u"Hz")

    tbl = uvtable(uv; stokes=UniversalSet)
    @test length(tbl) == 160560
    tbl = uvtable(uv)
    @test length(tbl) == 80280
    @test isconcretetype(eltype(tbl))
    @test tbl[12345] == (
        datetime = DateTime("2010-12-24T08:06:25"),
        stokes = :RR,
        freq_spec = VLBI.FrequencyWindow(6, 1.5369459f10u"Hz", 8.0f6u"Hz", 1, 1),
        spec = VLBI.VisSpec(VLBI.Baseline((:FD, :PT)), UV([4.282665f6, -2.8318948f7])),
        value = (0.48758677f0 - 0.09014242f0im) ±ᵤ 0.040410895f0
    )
    @test VLBI.antenna_names(tbl[12345]) == (:FD, :PT)
    @test VLBI.Baseline(tbl[12345]) == VLBI.Baseline((:FD, :PT))
    @test UV(tbl[12345]) == UV(4.282665f6, -2.8318948f7)
    @test VLBI.visibility(tbl[12345]) == (0.48758677f0 - 0.09014242f0im) ±ᵤ 0.040410895f0

    @test first(tbl) == first(uvtable(uv; impl=pyimport))
    @test tbl == uvtable(uv; impl=pyimport)
end

@testitem "closures calculations" begin
    # using VLBIDataExtra: read_uvtbl, closures_a_single, closures_a_all, closures_p_single, closures_p_all
    # import VLBIDataExtra: VLBIData
    using DataManipulation
    using VLBIFiles.Uncertain

    file = joinpath(dirname(pathof(VLBIFiles)), "../test/data/vis.fits")
    uvtbl = VLBI.load(file) |> uvtable

    # compute closures from data:
    gr = @p uvtbl group_vg((;_.freq_spec, _.stokes, _.datetime)) first
    @test length(gr) == 10

    cla_gr = VLBI.VLBIData.closures_scan(ClosureAmpSpec, gr)
    @test length(cla_gr) == 15
    @test cla_gr[7][[:stokes, :value,]] == (
        stokes = :RR,
        value = (0.7357933f0 + 0.08911123f0im) ±ᵤ 0.15352207f0,
    )

    clp_gr = VLBI.VLBIData.closures_scan(ClosurePhaseSpec, gr)
    @test length(clp_gr) == 10
    @test clp_gr[7][[:stokes, :value,]] == (
        stokes = :RR,
        value = (0.1127984f0 - 0.023669038f0im) ±ᵤ 0.021634521f0,
    )

    cla_all = VLBI.VLBIData.closures_all(ClosureAmpSpec, uvtbl)
    @test length(cla_all) == 507852
    @test (@p cla_all group_vg((;_.freq_spec, _.stokes, _.datetime)) first) == cla_gr

    clp_all = VLBI.VLBIData.closures_all(ClosurePhaseSpec, uvtbl)
    @test length(clp_all) == 146772
    @test (@p clp_all group_vg((;_.freq_spec, _.stokes, _.datetime)) first) == clp_gr
end


@testitem "uvf multichannel" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIFiles.Uncertain
    using Statistics
    using StaticArrays
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/vis_multichan.vis")
    @test length(uv.freq_windows) == 8
    df = uvtable(uv)
    target = (
        datetime = DateTime("1996-06-05T19:16:45.001"),
        stokes = :LL,
        freq_spec = VLBI.FrequencyWindow(3, 4.84099f9u"Hz", 8.0f6u"Hz", 16, 1),
        spec = VLBI.VisSpec(VLBI.Baseline((:HN, :LA)), UV([3.6239796f7, 1.2843347f7])),
        value = (-0.21484283f0 - 0.35979474f0im) ±ᵤ (1/√3.0233376f0)
    )
    res = filter(r -> VLBI.Baseline(r) == VLBI.Baseline((:HN, :LA)) && r.datetime == target.datetime && r.freq_spec == target.freq_spec && r.stokes == target.stokes, df)[1]
    @test res == target
    @test map(typeof, res) == map(typeof, target)
    @test df == uvtable(uv; impl=pyimport)
end

@testitem "uvf EHT 1" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIFiles.Uncertain
    using Statistics
    using StaticArrays
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/hops_3600_OJ287_LO+HI.medcal_dcal_full.uvfits")
    @test length(uv.freq_windows) == 2
    df = uvtable(uv)
    target = (
        datetime = Dates.DateTime("2017-04-10T00:58:45.005"), 
        stokes = :RR, 
        freq_spec = VLBI.FrequencyWindow(1, 2.270707f11u"Hz", 1.856f9u"Hz", 1, 1), 
        spec = VLBI.VisSpec(VLBI.Baseline((:AP, :AZ)), UV([2.9769375f9, -4.5123123f9])), 
        value = (0.5662228f0 + 0.014944182f0im) ±ᵤ (1/√117.818825f0)
    )
    res = filter(r -> VLBI.Baseline(r) == VLBI.Baseline((:AP, :AZ)) && r.datetime == target.datetime && r.freq_spec == target.freq_spec && r.stokes == target.stokes, df)[1]
    @test res == target
    @test map(typeof, res) == map(typeof, target)
    @test df == uvtable(uv; impl=pyimport)
end

@testitem "uvf EHT 2" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIFiles.Uncertain
    using Statistics
    using StaticArrays
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/SR1_3C279_2017_101_hi_hops_netcal_StokesI.uvfits")
    @test length(uv.freq_windows) == 1
    df = uvtable(uv)
    target = (
        datetime = Dates.DateTime("2017-04-11T02:14:55"), 
        stokes = :RR, 
        freq_spec = VLBI.FrequencyWindow(1, 2.290707f11u"Hz", 1.856f9u"Hz", 1, 1), 
        spec = VLBI.VisSpec(VLBI.Baseline((:AA, :AP)), UV([687115.25f0, -1.8692805f6])), 
        value = (-1.1796279f0 - 7.8919725f0im) ±ᵤ (1/√40441.19f0)
    )
    res = filter(r -> VLBI.Baseline(r) == VLBI.Baseline((:AA, :AP)) && r.datetime == target.datetime && r.freq_spec == target.freq_spec && r.stokes == target.stokes, df)[1]
    @test res == target
    @test map(typeof, res) == map(typeof, target)
    @test df == uvtable(uv; impl=pyimport)
end

@testitem "uvf EHT 3" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using Statistics
    using StaticArrays
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/datafile_01-01_230GHz.uvfits")
    @test length(uv.freq_windows) == 1
    df = uvtable(uv)
    @test VLBI.visibility(df[1]) ≈ -0.0061733737f0 + 0.052109245f0im
    @test df == uvtable(uv; impl=pyimport)
end

@testitem "difmap model" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    cd(dirname(@__FILE__))

    mod = VLBI.load("./data/difmap_model.mod")
    @test length(components(mod)) == 4
    @test map(flux, components(mod)) |> collect ≈ [0.522, 0.0217, 0.0176, 0.2145]u"Jy"  rtol=0.01
    @test map(fwhm_max, components(mod)) |> collect ≈ [0.135, 0.463, 1.99, 0]u"mas"  rtol=0.01
    @test coords(components(mod)[1]) ≈ [-0.000662, -0.00123]u"mas" rtol=0.01

    mod_empty = VLBI.load("./data/difmap_model_empty.mod")
    @test isempty(components(mod_empty))

    mod_clean = VLBI.load("./data/difmap_model_clean.mod")
    @test length(components(mod_clean)) == 631
    @test isconcretetype(eltype(components(mod_clean)))
    @test eltype(components(mod_clean)) <: Point

    mod_map = VLBI.load(MultiComponentModel, "./data/map.fits")

    tmpf = tempname()
    VLBI.save(tmpf, mod)
    @test VLBI.load(tmpf) == mod
    VLBI.save(tmpf, mod_empty)
    @test VLBI.load(tmpf) == mod_empty
    VLBI.save(tmpf, mod_clean)
    @test VLBI.load(tmpf) == mod_clean
    VLBI.save(tmpf, mod_map)
    @test VLBI.load(tmpf) ≈ mod_map  rtol=1e-4  # approx because saving to *.mod involves rounding; also Float32 vs 64
end

# @testitem "RFC" begin
#     using AstroRFC: RFC
#     using ProgressMeter
#     using PyCall

#     @testset "vis" begin
#         rfci = RFC.Files()
#         @showprogress for f in RFC.files(rfci, suffix="vis")
#             try
#                 uv = VLBI.load(f)
#                 df = uvtable(uv)
#                 @test df == uvtable(uv; impl=pyimport)
#             catch e
#                 @show abspath(f) e
#                 rethrow()
#             end
#         end
#     end

#     @testset "maps" begin
#         rfci = RFC.Files()
#         @showprogress for f in RFC.files(rfci, suffix="map", extension="fits")
#             try
#                 VLBI.load(f)
#             catch e
#                 @show f e
#                 rethrow()
#             end
#             try
#                 VLBI.load(MultiComponentModel, f)
#             catch e
#                 @show f e
#                 rethrow()
#             end
#             try
#                 VLBI.load(Beam, abspath(f))
#             catch e
#                 if e isa KeyError && e.key == "BMAJ"
#                     @warn "Missing BMAJ" f
#                 else
#                     @show f e
#                     rethrow()
#                 end
#             end
#         end
#     end
# end

@testitem "_" begin
    import Aqua
    Aqua.test_all(VLBIFiles; ambiguities=false, piracies=false)
    Aqua.test_ambiguities(VLBIFiles)

    import CompatHelperLocal as CHL
    CHL.@check()
end
