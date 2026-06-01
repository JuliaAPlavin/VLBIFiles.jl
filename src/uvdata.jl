"""
    uvw_wavelengths(uvw, freq) -> dimensionless wavelengths

Convert a UVW coordinate to wavelengths at the given frequency.
- `uvw::Unitful.Time` (FITS-IDI convention, light-travel time): returns `uvw * freq`.
- `uvw::Unitful.Length` (e.g. meters): returns `uvw / (c / freq)`.

Broadcast for vector quantities (e.g., a `UVW`).
"""
uvw_wavelengths(uvw::Unitful.Time,   freq::Unitful.Frequency) = NoUnits(uvw * freq)
# Use uvw / (c/freq) rather than uvw*freq/c — keeps the intermediate at moderate
# magnitude, avoiding Float32 precision loss for typical VLBI baselines.
uvw_wavelengths(uvw::Unitful.Length, freq::Unitful.Frequency) = NoUnits(uvw / (u"c" / freq))


Base.@kwdef struct FrequencyWindow
    freqid::Int
    ix::Int
    freq::typeof(1f0u"Hz")
    width::typeof(1f0u"Hz")
    nchan::Int16
    sideband::Int8
    crpix::Float32 = 1f0
end

frequency(fw::FrequencyWindow, kind::Symbol=:reference) = if kind == :reference
    fw.freq
elseif kind == :average
    @assert fw.nchan > 0
    cdelt = fw.sideband * fw.width / fw.nchan
    fw.freq + cdelt * ((fw.nchan + 1) / 2f0 - fw.crpix)
else
    error("Unsupported kind: $kind")
end

function frequency(fw::FrequencyWindow, ::Type{Interval})
    @assert fw.nchan > 0
    cdelt = fw.sideband * fw.width / fw.nchan
    lo = fw.freq + cdelt * (1 - fw.crpix) - cdelt / 2
    hi = fw.freq + cdelt * (fw.nchan - fw.crpix) + cdelt / 2
    return min(lo, hi) .. max(lo, hi)
end

frequency(freq::AbstractVector{<:VLBIFiles.FrequencyWindow}, ::Type{Interval}) = @p freq flatmap(endpoints(frequency(_, Interval))) extrema Interval(__...)

function frequencies(fw::FrequencyWindow)
    @assert fw.nchan > 0
    cdelt = fw.sideband * fw.width / fw.nchan
    return Float64(fw.freq) .+ ((1:fw.nchan) .- fw.crpix) .* cdelt
end

frequencies(fws::AbstractVector{<:FrequencyWindow}) = flatmap(frequencies, fws)

Base.isless(a::FrequencyWindow, b::FrequencyWindow) = a.freq < b.freq

function Statistics.mean(xs::AbstractVector{<:FrequencyWindow})
    nchan = (@p xs map(_.nchan) sum)
    FrequencyWindow(
        first(xs).freqid,
        first(xs).ix,
        (@p xs map(_.freq) mean),
        (@p xs map(_.width) sum),
        nchan,
        (@p xs map(_.sideband) uniqueonly),
        Float32((nchan + 1) / 2),
    )
end

function VLBIData._aggfreq(freq_specs::AbstractVector{<:FrequencyWindow})
    mean_freq = mean(fs -> fs.freq, freq_specs)
    total_width = sum(fs -> fs.width, freq_specs)
    return VLBIFiles.FrequencyWindow(0, 0, mean_freq, total_width, 1, 1, 1f0)
end

Base.@kwdef struct UVHeader
    fits::FITSHeader
    object::Union{String,Nothing}
    date_obs::Date
    stokes::Vector{Symbol}
    frequency::typeof(1.0u"Hz")
end

frequency(h::UVHeader) = h.frequency
Dates.Date(h::UVHeader) = h.date_obs
ICRSCoords(h::UVHeader) = ICRSCoords(h.fits["OBSRA"] * u"°", h.fits["OBSDEC"] * u"°")

function UVHeader(fh::FITSHeader)
    if (@oget fh["XTENSION"]) == "BINTABLE" && (@oget fh["EXTNAME"]) ∈ ["UV_DATA", "AIPS UV"]
        # FITS IDI
    else
        # uvfits
        axtypes = axis_types(fh)
        @assert "RA" in axtypes && "DEC" in axtypes
        @assert "FREQ" in axtypes
        @assert "COMPLEX" in axtypes && fh["NAXIS$(axis_ind(fh, "COMPLEX"))"] == 3
    end

    stokes = [val_to_stokes[val] for val in axis_vals(fh, "STOKES")]
    date_obs_str = fh["DATE-OBS"]
    date_obs = if occursin("/", date_obs_str)
        # Old FITS format: DD/MM/YY (20th century)
        d = Date(date_obs_str, dateformat"d/m/y")
        year(d) < 100 ? d + Year(1900) : d
    else
        # Modern format: YYYY-MM-DD (possibly with time suffix)
        Date(match(r"([\d-]+)", date_obs_str).captures[1], dateformat"Y-m-d")
    end

    return UVHeader(
        fits=fh,
        object=(@oget fh["OBJECT"]),
        date_obs=date_obs,
        stokes=stokes,
        frequency=axis_dict(fh, "FREQ")["CRVAL"]*u"Hz",
    )
end


function VLBIData.Antenna(hdu_row::NamedTuple)
    if haskey(hdu_row, :ORBPARM) && !isempty(hdu_row.ORBPARM) && hdu_row.ORBPARM != 0
        @warn "Antennas with ORBPARM detected, be careful" hdu_row.ORBPARM hdu_row.ANNAME
    end
    poltypes = haskey(hdu_row, :POLTYA) ? Symbol.((hdu_row.POLTYA, hdu_row.POLTYB)) : (:UNK, :UNK)
    # POLAA/POLAB are E(nband) per FITS-IDI v3; AIPS AN may carry scalars.
    # Use `first` to take the band-1 entry uniformly.
    feed_offsets = haskey(hdu_row, :POLAA) ?
        deg2rad.((Float64(first(hdu_row.POLAA)), Float64(first(hdu_row.POLAB)))) :
        (NaN, NaN)
    Antenna(; name=Symbol(hdu_row.ANNAME), xyz=hdu_row.STABXYZ, mount_type=VLBIData.AntennaMountType.T(hdu_row.MNTSTA), poltypes, feed_offsets)
end
Base.@kwdef struct AntArray
    name::String
    freq::typeof(1f0u"Hz")
    antennas::Dictionary{Int, Antenna}
    array_xyz::SVector{3, Float64} = SVector(NaN, NaN, NaN)
end

strfloat_to_float(x::Real) = x
strfloat_to_float(x::String) = parse(Float64, replace(x, "D" => "E"))

function AntArray(hdu::TableHDU)
    header = read_header(hdu)
    antennas = map(row -> row.NOSTA => Antenna(row), hdu |> columntable |> StructArray) |> dictionary
    array_xyz = SVector(strfloat_to_float.((
        get(header, "ARRAYX", NaN),
        get(header, "ARRAYY", NaN),
        get(header, "ARRAYZ", NaN),
    )))
    AntArray(;
        name=header["ARRNAM"],
        freq=strfloat_to_float(header["FREQ"]) * u"Hz",
        antennas,
        array_xyz,
    )
end

Base.length(a::AntArray) = length(a.antennas)
Base.getindex(a::AntArray, i::Int) = a.antennas[i]
Base.setindex(a::AntArray, ant::Antenna, i::Int) = @set a.antennas[i] = ant
ant_by_name(a::AntArray, name::Symbol) = filteronly(ant -> ant.name == name, a.antennas)


function Baseline_from_fits(b::Real)
    bi = floor(Int, b)
    a1, a2 = bi ÷ 256, bi % 256
    return Baseline((a1, a2))
end
function Baseline_from_fits(b::Real, antarrays)
    bi = floor(Int, b)
    arri = round(Int, (b % 1) * 100) + 1
    a1, a2 = bi ÷ 256, bi % 256
    return Baseline(arri, (a1, a2), antarrays)
end

function VLBIData.Baseline(array_ix::Integer, ant_ids::NTuple{2, Integer}, ant_arrays::Vector{AntArray})
    ants = ant_arrays[array_ix].antennas
    names = map(ant_ids) do id
        ant = get(ants, id, nothing)
        if !isnothing(ant)
            ant.name
        else
            @warn "Antenna index out of bounds, assigning generated name" length(ants) ant_ids maxlog=1
            Symbol(:ANT, id)
        end
    end
    Baseline(names)
end


DateTime_from_DATE_TIME(DATE, TIME) = julian_day(Float64(DATE) + TIME)


Base.@kwdef struct UVData
    path::String
    header::UVHeader
    freq_windows::Vector{FrequencyWindow}
    ant_arrays::Vector{AntArray}
end

FITSIO.FITS(uvdata::UVData) = FITS(uvdata.path)

function read_freqs(uvh, fq_table; crpix)
    ref_freq = @oget frequency(uvh) read_header(fq_table)["REF_FREQ"]*u"Hz"
	fq_table = fq_table |> StructArrays.fromtable
    res = @p fq_table flatmap() do fq_row
		subrows = if any(x->x isa AbstractVector, fq_row)
			@p fq_row |> filter(_ isa AbstractVector) |> rowtable
		else
			[fq_row]
		end
		@p subrows |> enumerate |> map() do (ix, r)
	        total_bw = @oget r[S"TOTAL BANDWIDTH"] r[S"TOTAL_BANDWIDTH"]
	        ch_width = @oget r[S"CH WIDTH"] r[S"CH_WIDTH"]
	        curfreq = @oget r[S"IF FREQ"] r[S"BANDFREQ"]
	        nchan = Int(total_bw / abs(ch_width))
	        FrequencyWindow(;
                freqid=(@oget fq_row.FREQID 1),
	            ix,
	            freq=ref_freq + curfreq * u"Hz",
	            width=total_bw * u"Hz",
	            nchan,
	            sideband=r.SIDEBAND,
	            crpix,
	        )
	    end
	end
end


function uvw_keys(colnames)
    tag(c) = lstrip(string(c)[3:end], '-')          # suffix after the 2-letter prefix and any dashes
    uu = only(filter(c -> startswith(string(c), "UU"), colnames))
    t  = tag(uu)
    pick(pre) = only(filter(c -> startswith(string(c), pre) && tag(c) == t, colnames))
    (uu, pick("VV"), pick("WW"))
end


function to_proper_stokes(stokes::Symbol, ant_arrays, bl::Baseline{Symbol})
    length(ant_arrays) == 1 || error("Modify poltypes in multi-array UVData is not supported")
    ant_array = only(ant_arrays)
    ants = ant_by_name.(Ref(ant_array), bl.antennas)
    poltypes_is =
        stokes == :RR ? (1, 1) :
        stokes == :LL ? (2, 2) :
        stokes == :RL ? (1, 2) :
        stokes == :LR ? (2, 1) :
        error("Unsupported stokes: $stokes")
    poltypes = map((a, i) -> a.poltypes[i], ants, poltypes_is)
    return Symbol(poltypes...)
end


function uvtable(uv::UVData)
    wt = uvtable_wide(uv)
    stk = uv.header.stokes
    nch = uniqueonly(fw.nchan for fw in uv.freq_windows)
    # stokes_asis: when every antenna shares the same feeds, the file's stokes labels already match the
    # data, so skip the per-baseline feed translation; uniform-feed files (including linear/Stokes) use
    # their raw labels, and only mixed-feed arrays are translated.
    stokes_asis = length(@p uv.ant_arrays flatmap(_.antennas) map(_.poltypes) unique) == 1
    properstokes(s, bl) = stokes_asis ? s : to_proper_stokes(s, uv.ant_arrays, bl)
    @p begin
        wt
        flatmap() do r
            # one row per surviving (stokes, channel); `product` varies its first argument fastest, so the
            # order is IF → channel-within-IF → stokes-innermost. weight ≤ 0 is dropped; k and w computed once.
            @p Iterators.product(eachindex(stk), 1:nch, enumerate(uv.freq_windows)) filtermap() do (s, chan, (iif, fw))
                k = chan + nch*(iif - 1)
                w = r.weight[s, k]
                w > 0 || return nothing
                (; r.datetime,
                   stokes = properstokes(stk[s], r.baseline),
                   freq_spec = fw,
                   spec = VisSpec(r.baseline, UV(UVW(uvw_wavelengths.(r.uvw, frequency(fw, :average))...))),
                   value = U.Value(r.visibility[s, k], 1/√w))
            end
        end
        StructArray()
    end
end


function uvtable_wide(uv::UVData)
    FITS(uv.path) do fits
        haskey(fits, "UV_DATA") ? _widetable_fitsidi(uv, fits) : _widetable_uvfits(uv, fits)
    end
end

ICRSCoords(uv::UVData) = ICRSCoords(uv.header)

function sources(uv::UVData)
    FITS(uv.path) do fits
        haskey(fits, "SOURCE") ? _sources_fitsidi(fits) : _sources_uvfits(uv)
    end
end

_sources_uvfits(uv::UVData) = dictionary([1 => (; name = uv.header.object, coords = ICRSCoords(uv.header))])

function _sources_fitsidi(fits)
    rows = StructArrays.fromtable(fits["SOURCE"])
    dictionary(Int(@oget r.SOURCE_ID r.var"ID_NO.") =>
               (name = String(strip(r.SOURCE)), coords = ICRSCoords(r.RAEPO * u"°", r.DECEPO * u"°"))
               for r in rows)
end

# shared per-row helpers; `cols` is a column table (lazy or materialized).
_decode_uvw(cols)          = (k = uvw_keys(keys(cols)); map((u,v,w) -> SVector(u,v,w) .* u"s" .* u"c" .|> u"m",
                                                            cols[k[1]], cols[k[2]], cols[k[3]]))
_layout(uv)               = (length(uv.header.stokes),                        # n_stokes
                             uniqueonly(fw.nchan for fw in uv.freq_windows),   # nchan_per_IF
                             length(uv.freq_windows))                          # n_IF

# FITS-IDI UVW are light-travel seconds (TUNIT may be absent in older files); refuse any other unit.
function _assert_uvw_seconds(fh, ks)
    units = Dict(strip(fh["TTYPE$i"]) => (@oget fh["TUNIT$i"]) for i in 1:fh["TFIELDS"])
    for k in ks
        u = units[string(k)]
        isnothing(u) || strip(u) ∈ ("SECONDS", "s", "sec") || error("FITS-IDI UVW column $k has TUNIT=$u, expected seconds")
    end
end

# L1 analysis cell: transform one record's L0 N-D array `nd` into a (stokes,channel) complex KeyedArray.
# `keyless_unname` hands FitsVMatrix the raw block so it can reshape COMPLEX/stokes/(FREQ×BAND) into its
# named (COMPLEX,stokes,channel) view (re/im combined on access). The stokes axis carries over from `nd`
# (same values); only `freq` is supplied fresh — the authoritative per-channel axis from freq_windows.
_viscell(nd, nc, ns, nt, frq) =
    KeyedArray(FitsVMatrix(keyless_unname(nd), nc, ns, nt); stokes=axiskeys(nd, :STOKES), freq=frq)

function _widetable_fitsidi(uv, fits)
    length(unique(fw.freqid for fw in uv.freq_windows)) == 1 || error("multi-FREQID files not supported")
    hdu  = fits["UV_DATA"]
    fh   = read_header(hdu)
    cols = lazycolumntable(hdu)                             # L0: FLUX/WEIGHT are faithful N-D KeyedArray columns
    _assert_uvw_seconds(fh, uvw_keys(keys(cols)))
    n_stokes, nchan_per_IF, n_IF = _layout(uv); n_total = nchan_per_IF * n_IF
    n_complex = fh["MAXIS1"]
    # eager scalar metadata in one pread pass (FLUX/WEIGHT — the array columns — are dropped first)
    meta = materialize_columns(delete(cols, @o _.FLUX _.WEIGHT))
    baseline  = map(b -> Baseline_from_fits(b, uv.ant_arrays), meta.BASELINE)
    datetime  = julian_day.(meta.DATE .+ meta.TIME)
    uvw       = _decode_uvw(meta)
    source_id = Int.(meta.SOURCE)
    frq = frequencies(uv.freq_windows)                      # authoritative per-channel freq; stokes reused from each cell
    # L1: each cell is a `mapview` transform of the L0 N-D array — no `parent`, the lazy column is `mapview` itself
    visibility = mapview(nd -> _viscell(nd, n_complex, n_stokes, n_total, frq), cols.FLUX)
    weight = if haskey(cols, :WEIGHT)                        # separate WEIGHT N-D (per stokes, band)
        mapview(w -> KeyedArray(FitsWMatrixSep(keyless_unname(w), n_stokes, nchan_per_IF, n_total); stokes=axiskeys(w, :STOKES), freq=frq), cols.WEIGHT)
    else                                                     # embedded weight = 3rd COMPLEX entry of FLUX
        mapview(nd -> KeyedArray(FitsWMatrixEmb(keyless_unname(nd), n_complex, n_stokes, n_total); stokes=axiskeys(nd, :STOKES), freq=frq), cols.FLUX)
    end
    wt = StructArray((; baseline, datetime, uvw, source_id, visibility, weight))
    _finalize(wt, meta, uv)
end

function _widetable_uvfits(uv, fits)
    hdu  = GroupedHDU(fits.fitsfile, 1)
    cols = lazycolumntable(hdu)                             # L0: DATA is a faithful N-D KeyedArray column
    isnothing(cols) && error("UVFITS file is not mmappable")
    n_stokes, nchan_per_IF, n_IF = _layout(uv); n_total = nchan_per_IF * n_IF
    n_complex = read_header(hdu)["NAXIS2"]
    n_complex == 3 || error("UVFITS COMPLEX axis NAXIS2=$n_complex; only embedded-weight files (NAXIS2=3) are supported")
    # eager scalar metadata (the DATA array column is dropped first, stays lazy) — same call as FITS-IDI
    meta = materialize_columns(delete(cols, @o _.DATA))
    baseline  = map(b -> Baseline_from_fits(b, uv.ant_arrays), meta.BASELINE)
    date2     = haskey(meta, :_DATE) ? meta._DATE : 0
    # the two-word JD (DATE integer part + _DATE fraction) must be summed in Float64: both words
    # are Float32 on disk, but their sum (~2.45e6 days) has no Float32 representation.
    datetime  = julian_day.(Float64.(meta.DATE) .+ date2)
    uvw       = _decode_uvw(meta)
    source_id = fill(1, length(baseline))
    frq = frequencies(uv.freq_windows)
    visibility = mapview(nd -> _viscell(nd, n_complex, n_stokes, n_total, frq), cols.DATA)
    weight     = mapview(nd -> KeyedArray(FitsWMatrixEmb(keyless_unname(nd), n_complex, n_stokes, n_total); stokes=axiskeys(nd, :STOKES), freq=frq), cols.DATA)
    wt = StructArray((; baseline, datetime, uvw, source_id, visibility, weight))
    _finalize(wt, meta, uv)
end

# adds the optional int_time column and emits faithful warnings; `cols` is the row-metadata column table
function _finalize(wt, cols, uv)
    haskey(cols, :INTTIM) && (wt = @insert wt.int_time = cols.INTTIM .* u"s")
    issorted(wt.datetime)                   || @warn "visibility rows are not time-sorted" path=uv.path
    _has_noncanonical_baseline(wt.baseline) && @warn "some physical baselines appear in both orientations" path=uv.path
    wt
end

# §8 "canonical" = each physical baseline uses one orientation throughout (FringyFITS asserts this); flag
# if some {a,b} appears as both (a,b) and (b,a). Reads the already-decoded baselines (no raw-code re-decode).
function _has_noncanonical_baseline(baselines)
    prs = Set(VLBIData.antenna_names(bl) for bl in baselines)
    any(p -> p[1] != p[2] && (p[2], p[1]) ∈ prs, prs)
end

function load(::Type{UVData}, path)
    path = abspath(path)  # for RFC.File
    fits = FITS(path)
    header = let
        fh = if haskey(fits, "UV_DATA")
            read_header(fits["UV_DATA"])  # fits idi
        elseif haskey(fits, "AIPS UV")
            read_header(fits["AIPS UV"])  # VLA uvfits?.. don't really work
        else
            read_header(fits[1])  # uvfits
        end
        UVHeader(fh)
    end
    freq_crpix = get(axis_dict(header.fits, "FREQ"), "CRPIX", 1.0)
    fq_table = @oget fits["AIPS FQ"] fits["FREQUENCY"] nothing
    freq_windows = if fq_table !== nothing
        read_freqs(header, fq_table; crpix=freq_crpix)
    else
        # No FQ table: construct from header FREQ axis
        fd = axis_dict(header.fits, "FREQ")
        nchan = @oget fd["NAXIS"] fd["MAXIS"]
        cdelt = fd["CDELT"]
        [FrequencyWindow(;
            freqid=1,
            ix=1,
            freq=Float32(fd["CRVAL"]) * u"Hz",
            width=Float32(nchan * abs(cdelt)) * u"Hz",
            nchan=Int16(nchan),
            sideband=Int8(sign(cdelt)),
            crpix=Float32(freq_crpix),
        )]
    end
    ant_arrays = AntArray[]
    for i in Iterators.countfrom(1)
        hdu = try
            haskey(fits, "AIPS AN") ? fits["AIPS AN", i] : fits["ARRAY_GEOMETRY", i]
        catch err
            err isa FITSIO.CFITSIO.CFITSIOError && "illegal HDU number" == err.errmsgshort && break
            rethrow()
        end
        push!(ant_arrays, AntArray(hdu))
    end

    # Per FITS-IDI v3 (AIPS Memo 114r), the global ANTENNA HDU carries per-antenna
    # feed polarization (POLTYA/POLTYB) and feed-offset angles (POLAA/POLAB).
    # ARRAY_GEOMETRY does not, so we merge here.
    # POLAA/POLAB are E(nband) per FITS-IDI v3 (per-row vector of length nband);
    # AIPS-style scalars come through as a plain number. `first` picks the band-1
    # entry uniformly in either case.
    if haskey(fits, "ANTENNA")
        at_rows = fits["ANTENNA"] |> columntable |> StructArray
        for ai in eachindex(ant_arrays)
            arr = ant_arrays[ai]
            for (id, ant) in pairs(arr.antennas)
                idx = findfirst(r -> r.ARRAY == ai && r.ANTENNA_NO == id, at_rows)
                idx === nothing && continue
                r = at_rows[idx]
                arr.antennas[id] = Antenna(;
                    name         = ant.name,
                    xyz          = ant.xyz,
                    mount_type   = ant.mount_type,
                    poltypes     = (Symbol(strip(r.POLTYA)), Symbol(strip(r.POLTYB))),
                    feed_offsets = deg2rad.((Float64(first(r.POLAA)),
                                              Float64(first(r.POLAB)))),
                )
            end
        end
    end

    close(fits)

    UVData(; path, header, freq_windows, ant_arrays)
end
