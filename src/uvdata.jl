Base.@kwdef struct FrequencyWindow
    ix::Int
    freq::typeof(1f0u"Hz")
    width::typeof(1f0u"Hz")
    nchan::Int16
    sideband::Int8
end

frequency(fw::FrequencyWindow, kind::Symbol=:reference) = if kind == :reference
    fw.freq
elseif kind == :average
    @assert fw.sideband == 1
    @assert fw.nchan > 0
    fw.freq + fw.width / 2
else
    error("Unsupported kind: $kind")
end

function frequency(fw::FrequencyWindow, ::Type{Interval})
    @assert fw.sideband == 1
    @assert fw.nchan > 0
    return fw.freq .. (fw.freq + fw.width)
end

function frequencies(fw::FrequencyWindow)
    @assert fw.sideband == 1
    @assert fw.nchan > 0
    return Float64(fw.freq) .+ (0.5:(fw.nchan-0.5)) .* (fw.width/fw.nchan)
end

Base.isless(a::FrequencyWindow, b::FrequencyWindow) = a.freq < b.freq

Statistics.mean(xs::AbstractVector{<:FrequencyWindow}) = FrequencyWindow(
    first(xs).ix,
    (@p xs map(_.freq) mean),
    (@p xs map(_.width) sum),
    (@p xs map(_.nchan) sum),
    (@p xs map(_.sideband) uniqueonly),
)

Base.@kwdef struct UVHeader
    fits::FITSHeader
    object::Union{String,Nothing}
    date_obs::Date
    stokes::Vector{Symbol}
    frequency::typeof(1.0u"Hz")
end

frequency(h::UVHeader) = h.frequency
Dates.Date(h::UVHeader) = h.date_obs

function UVHeader(fh::FITSHeader)
    if (@oget fh["XTENSION"]) == "BINTABLE" && (@oget fh["EXTNAME"]) == "UV_DATA"
        # FITS IDI
    else
        # uvfits
        @assert fh["CTYPE6"] == "RA" && fh["CTYPE7"] == "DEC"
        @assert fh["CTYPE4"] == "FREQ"
        @assert fh["CTYPE2"] == "COMPLEX" && fh["NAXIS2"] == 3
    end

    stokes = [val_to_stokes[val] for val in axis_vals(fh, "STOKES")]
    date_obs = match(r"([\d-]+)(\(\d+\))?", fh["DATE-OBS"]).captures[1]

    return UVHeader(
        fits=fh,
        object=(@oget fh["OBJECT"]),
        date_obs=Date(date_obs, dateformat"Y-m-d"),
        stokes=stokes,
        frequency=axis_dict(fh, "FREQ")["CRVAL"]*u"Hz",
    )
end


function VLBIData.Antenna(hdu_row::NamedTuple)
    if !isempty(hdu_row.ORBPARM) && hdu_row.ORBPARM != 0
        @warn "Antennas with ORBPARM detected, be careful" hdu_row.ORBPARM hdu_row.ANNAME
    end
    poltypes = haskey(hdu_row, :POLTYA) ? Symbol.((hdu_row.POLTYA, hdu_row.POLTYB)) : (:UNK, :UNK)
    Antenna(; name=Symbol(hdu_row.ANNAME), xyz=hdu_row.STABXYZ, mount_type=VLBIData.AntennaMountType.T(hdu_row.MNTSTA), poltypes)
end
Base.@kwdef struct AntArray
    name::String
    freq::typeof(1f0u"Hz")
    antennas::Dictionary{Int, Antenna}
end

strfloat_to_float(x::AbstractFloat) = x
strfloat_to_float(x::String) = parse(Float64, replace(x, "D" => "E"))

function AntArray(hdu::TableHDU)
    header = read_header(hdu)
    antennas = map(row -> row.NOSTA => Antenna(row), hdu |> columntable |> StructArray) |> dictionary
    AntArray(;
        name=header["ARRNAM"],
        freq=strfloat_to_float(header["FREQ"]) * u"Hz",
        antennas,
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
    arri = floor(Int, (b % 1) * 100) + 1
    a1, a2 = bi ÷ 256, bi % 256
    return Baseline(1, (a1, a2), antarrays)
end

function VLBIData.Baseline(array_ix::Integer, ant_ids::NTuple{2, Integer}, ant_arrays::Vector{AntArray})
    ants = ant_arrays[array_ix].antennas
    names = map(ant_ids) do id
        ant = get(ants, id, nothing)
        if !isnothing(ant)
            ant.name
        else
            @warn "Antenna index out of bounds, assigning generated name" length(ants) ant_ids
            Symbol(:ANT, id)
        end
    end
    Baseline(names)
end


Base.@kwdef struct UVData
    path::String
    header::Union{UVHeader,Nothing}
    freq_windows::Vector{FrequencyWindow}
    ant_arrays::Vector{AntArray}
end

FITSIO.FITS(uvdata::UVData) = FITS(uvdata.path)

function read_freqs(uvh, fq_table)
    fq_row = fq_table |> columntable |> StructArray |> only
    # fq_row = fq_row[Symbol.(["IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND"])]
    nrows = @p fq_row values() filter(_ isa AbstractVector) (isempty(__) ? 1 : length(__[1]))
    fq_row = map(fq_row) do x
        isa(x, AbstractArray) ? x : fill(x, nrows)
    end
    ref_freq = @oget frequency(uvh) read_header(fq_table)["REF_FREQ"]*u"Hz"
    res = map(fq_row |> rowtable |> enumerate) do (ix, r)
        total_bw = @oget r[S"TOTAL BANDWIDTH"] r[S"TOTAL_BANDWIDTH"]
        ch_width = @oget r[S"CH WIDTH"] r[S"CH_WIDTH"]
        curfreq = @oget r[S"IF FREQ"] r[S"BANDFREQ"]
        nchan = Int(total_bw / ch_width)
        FrequencyWindow(;
            ix,
            freq=ref_freq + curfreq * u"Hz",
            width=total_bw * u"Hz",
            nchan,
            sideband=r.SIDEBAND,
        )
    end
end


function read_data_raw(uvdata::UVData, ::typeof(identity)=identity)
    fits = FITS(uvdata)
    if haskey(fits, "UV_DATA")
        return fits["UV_DATA"] |> lazycolumntable |> StructArray
    end
    hdu = GroupedHDU(fits.fitsfile, 1)
    read(hdu)
end


function read_data_arrays(uvdata::UVData, impl=identity)
    raw = read_data_raw(uvdata, impl)

    uvw_keys = @p begin
        [(S"UU", S"VV", S"WW"), (S"UU--", S"VV--", S"WW--"), (S"UU---SIN", S"VV---SIN", S"WW---SIN")]
        filteronly(_ ⊆ keys(raw))
    end

    count = size(raw[:DATA])[end]
    n_if = length(uvdata.freq_windows)
    n_chan = map(fw -> fw.nchan, uvdata.freq_windows) |> unique |> only
    axarr = KeyedArray(raw[:DATA],
        COMPLEX=[:re, :im, :wt],
        STOKES=uvdata.header.stokes,
        FREQ=1:n_chan,
        IF=1:n_if,
        RA=[1], DEC=[1],
        _=1:count,
    )
    # drop always-singleton axes
    axarr = axarr[DEC=1, RA=1]

    uvw_m = UVW.([Float32.(raw[k]) .* (u"c" * u"s") .|> u"m" for k in uvw_keys]...)
    baseline = map(b -> Baseline_from_fits(b, uvdata.ant_arrays), raw[:BASELINE])
    data = (;
        uvw_m,
        baseline,
        datetime = julian_day.(Float64.(raw[:DATE]) .+ raw[:_DATE]),
        visibility = complex.(axarr(COMPLEX=:re), axarr(COMPLEX=:im)),
        weight = axarr(COMPLEX=:wt),
    )
    if haskey(raw, :INTTIM)
        data = merge(data, (int_time = raw[:INTTIM] .* u"s",))
    end
    return data
end

function _table(uvdata::UVData; impl)
    data = read_data_arrays(uvdata, impl)
    @assert ndims(data.visibility) == 4

    poltypes = @p uvdata.ant_arrays flatmap(_.antennas) map(_.poltypes) unique
    # XXX: will also trigger if poltypes of all antennas changed to a different (but same) value:
    stokes_asis = length(poltypes) == 1
    
    @p begin
        data.visibility
        # `|> columntable |> rowtable` is faster than `|> rowtable` alone
        map(__ |> columntable |> rowtable) do r
            ix = r._
            freq_spec = uvdata.freq_windows[r.IF]
            uvw_m = data.uvw_m[ix]
            uvw_wl = UVW(ustrip.(Unitful.NoUnits, uvw_m ./ (u"c" / frequency(freq_spec, :average))))
            (;
                baseline=data.baseline[ix],
                datetime=data.datetime[ix],
                stokes=stokes_asis ? r.STOKES : to_proper_stokes(r.STOKES, uvdata.ant_arrays, data.baseline[ix]),
                freq_spec,
                uv_m=UV(uvw_m), w_m=uvw_m.w,
                uv=UV(uvw_wl), w=uvw_wl.w,
                visibility=r.value,
            )
        end
        StructArray()
        @insert __.weight = collect(vec(data.weight))
        filter!(_.weight > 0)
    end
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


uvtable(uvd::UVData; impl=identity) = @p uvd _table(;impl) map((;
    _.datetime, _.stokes, _.freq_spec,
    spec=VisSpec(_.baseline, UV(_.uv)),
    value=U.Value(_.visibility, 1/√_.weight),
))

function load(::Type{UVData}, path)
    path = abspath(path)  # for RFC.File
    fits = FITS(path)
    header = let
        fh = read_header(fits[1])
        if haskey(fh, "CORRELAT")
            # fits idi
            fh = read_header(fits["UV_DATA"])
        end
        try
            UVHeader(fh)
        catch e
            e isa KeyError || rethrow()
            rethrow()
            nothing
        end
    end
    freq_windows = read_freqs(header, @oget fits["AIPS FQ"] fits["FREQUENCY"])
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
    close(fits)

    UVData(; path, header, freq_windows, ant_arrays)
end
