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

Statistics.mean(xs::AbstractVector{<:FrequencyWindow}) = FrequencyWindow(
    first(xs).ix,
	(@p xs map(_.freq) mean),
	(@p xs map(_.width) sum),
	(@p xs map(_.nchan) sum),
	(@p xs map(_.sideband) uniqueonly),
)

Base.@kwdef struct UVHeader
    fits::FITSHeader
    object::String
    date_obs::Date
    stokes::Vector{Symbol}
    frequency::typeof(1.0u"Hz")
end

frequency(h::UVHeader) = h.frequency
Dates.Date(h::UVHeader) = h.date_obs

function UVHeader(fh::FITSHeader)
    @assert fh["CTYPE6"] == "RA" && fh["CTYPE7"] == "DEC"
    @assert fh["CTYPE4"] == "FREQ"
    @assert fh["CTYPE2"] == "COMPLEX" && fh["NAXIS2"] == 3

    val_to_stokes = Dict(-4=>:LR, -3=>:RL, -2=>:LL, -1=>:RR, 1=>:I, 2=>:Q, 3=>:U, 4=>:V)
    stokes = [val_to_stokes[val] for val in axis_vals(fh, "STOKES")]
    date_obs = match(r"([\d-]+)(\(\d+\))?", fh["DATE-OBS"]).captures[1]

    return UVHeader(
        fits=fh,
        object=fh["OBJECT"],
        date_obs=Date(date_obs, dateformat"Y-m-d"),
        stokes=stokes,
        frequency=axis_dict(fh, "FREQ")["CRVAL"]*u"Hz",
    )
end


function VLBIData.Antenna(hdu_row::NamedTuple)
    if !isempty(hdu_row.ORBPARM) && hdu_row.ORBPARM != 0
        @warn "Antennas with ORBPARM detected, be careful" hdu_row.ORBPARM hdu_row.ANNAME
    end
    Antenna(; name=Symbol(hdu_row.ANNAME), xyz=hdu_row.STABXYZ)
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
    fits = FITS(uvdata.path)
    if haskey(fits, "UV_DATA")
        return (fits["UV_DATA"] |> columntable)
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
    baseline = map(raw[:BASELINE]) do b
        bi = floor(Int, b)
        Baseline(round(Int, (b % 1) * 100) + 1, (bi ÷ 256, bi % 256), uvdata.ant_arrays)
    end
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

function _table(uvdata::UVData, impl=identity)
    data = read_data_arrays(uvdata, impl)
    @assert ndims(data.visibility) == 4
    
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
                stokes=r.STOKES,
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


uvtable(uvd::UVData; stokes=(:I, :LL, :RR)) = @p uvd _table filter(_.stokes ∈ stokes) map((;
	_.datetime, _.stokes, _.freq_spec,
	spec=VisSpec(_.baseline, UV(_.uv)),
	value=U.Value(_.visibility, 1/√_.weight),
))

function load(::Type{UVData}, path)
    path = abspath(path)  # for RFC.File
    fits = FITS(path)
    fh = read_header(fits[1])
    header = try
        UVHeader(fh)
    catch e
        e isa KeyError || rethrow()
        nothing
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
