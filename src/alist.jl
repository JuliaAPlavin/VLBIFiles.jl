struct Alist end

function load(::Type{Alist}, file)
    raw = CSV.read(file, StructArrays.fromtable; ntasks=1, comment="*", delim=' ', ignorerepeated=true,
        header=["version", "root_id", "two", "extent_no", "duration", "length", "offset", "expt_no", "scan_id", "procdate", "year", "timetag", "scan_offset", "source", "baseline", "quality", "freq_code", "polarization", "lags", "amp", "snr", "resid_phas", "phase_snr", "datatype", "sbdelay", "mbdelay", "ambiguity", "delay_rate", "ref_elev", "rem_elev", "ref_az", "rem_az", "u", "v", "esdesp", "epoch", "ref_freq", "total_phas", "total_rate", "total_mbdelay", "total_sbresid", "srch_cotime", "noloss_cotime", "ra_hrs", "dec_deg", "resid_delay"],
        stringtype=String)
    @assert all(==(2), raw.two)
    map(raw) do x
        m = match(r"^(?<doy>\d\d\d)-(?<h>\d\d)(?<m>\d\d)(?<s>\d\d)$", x.timetag)
        datetime = Date(x.year) + Day(parse(Int, m[:doy]) - 1) + Time(parse(Int, m[:h]), parse(Int, m[:m]), parse(Int, m[:s]))
        freq_spec = x.ref_freq*u"MHz"
        (;
            x.root_id, x.source, datetime, x.scan_id, x.length,
            freq_spec, spec=VisSpec(Baseline(Symbol.(Tuple(x.baseline))), UV(x.u, x.v)),
            stokes=Symbol(x.polarization),
            value=x.amp / 1e4 * cis(deg2rad(x.resid_phas)) * (1 ±ᵤ (1 / x.snr)),
            x.snr
        )
    end
end
