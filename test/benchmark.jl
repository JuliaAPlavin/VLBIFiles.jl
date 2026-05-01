# julia --project=test test/benchmark.jl
using Random
using VLBIFiles

const DEFAULT_PATH = "/Users/aplavin/work/galactic_scatter/2005+403/data/archival/vlba/UG002/VLBA_UG002O_ug002o_BIN0_SRC0_0_180821T142741.idifits"
# const DEFAULT_PATH = joinpath(@__DIR__, "data/BL146_1.fits")
const SINK = Ref{Any}()

bench_path() = get(ENV, "VLBIFILES_BENCHMARK_PATH", DEFAULT_PATH)

function min_timed(f; samples=7)
    times = Float64[]
    bytes = Int[]
    SINK[] = f()
    for _ in 1:samples
        GC.gc()
        t = @timed (SINK[] = f())
        push!(times, t.time)
        push!(bytes, t.bytes)
    end
    (; time=minimum(times), bytes=minimum(bytes))
end

function report(label, nops, stats; nrows=nothing)
    us = stats.time / nops * 1e6
    bytes = stats.bytes / nops
    file_s = isnothing(nrows) ? "" : "  $(lpad(round(stats.time / nops * nrows, digits=2), 8)) s/file"
    println("  $(rpad(label, 34)) $(lpad(round(us, digits=3), 10)) us/op  $(lpad(round(bytes, digits=1), 10)) B/op$file_s")
end

function bench(label, nops, f; nrows=nothing, samples=7)
    report(label, nops, min_timed(f; samples); nrows)
end

function make_indices(nrows, n; rng=MersenneTwister(0x5151))
    n = min(n, nrows)
    random = rand(rng, 1:nrows, n)
    sorted = sort(copy(random))
    first_contiguous = rand(rng, 1:max(nrows - n + 1, 1))
    contiguous = collect(first_contiguous:(first_contiguous + n - 1))
    same_row = fill(rand(rng, 1:nrows), n)
    (; random, sorted, contiguous, same_row)
end

index_sets(indices) = (
    ("random", indices.random),
    ("sorted", indices.sorted),
    ("contiguous", indices.contiguous),
    ("same-row", indices.same_row),
)

function each_index(f, indices)
    result = nothing
    for i in indices
        result = f(i)
    end
    return result
end

function sampled_rows(raw, indices)
    rows = Vector{eltype(raw)}(undef, length(indices))
    for (j, i) in pairs(indices)
        rows[j] = raw[i]
    end
    return rows
end

p = bench_path()
uvf = VLBI.load(VLBI.UVData, p)

t = @elapsed raw = VLBI.read_data_raw(uvf)
nrows = length(raw)
println("file: $p")
println("read_data_raw: $(round(t, digits=2))s  ($nrows rows)")

n = parse(Int, get(ENV, "VLBIFILES_BENCHMARK_N", "10000"))
samples = parse(Int, get(ENV, "VLBIFILES_BENCHMARK_SAMPLES", "7"))
indices = make_indices(nrows, n)
println("elapsed/allocation minima over $samples samples")
println("  $(rpad("benchmark", 34)) $(lpad("elapsed", 16))  $(lpad("alloc", 15))")

println("\nscalar collect")
bench("collect SOURCE", nrows, () -> collect(raw.SOURCE); nrows=nrows, samples=samples)
bench("collect DATE", nrows, () -> collect(raw.DATE); nrows=nrows, samples=samples)

println("\nscalar access")
for (label, ix) in index_sets(indices)
    bench("SOURCE[$label]", length(ix), () -> each_index(i -> raw.SOURCE[i], ix); nrows=nrows, samples=samples)
end

println("\nFLUX row access")
for (label, ix) in index_sets(indices)
    bench("FLUX[$label]", length(ix), () -> each_index(i -> raw.FLUX[i], ix); nrows=nrows, samples=samples)
end

println("\nFLUX first element")
for (label, ix) in index_sets(indices)
    bench("FLUX[$label][1]", length(ix), () -> each_index(i -> raw.FLUX[i][1], ix); nrows=nrows, samples=samples)
end

println("\nsum(FLUX row)")
for (label, ix) in index_sets(indices)
    bench("sum(FLUX[$label])", length(ix), () -> each_index(i -> sum(raw.FLUX[i]), ix); nrows=nrows, samples=samples)
end

println("\nWEIGHT first element")
for (label, ix) in index_sets(indices)
    bench("WEIGHT[$label][1]", length(ix), () -> each_index(i -> raw.WEIGHT[i][1], ix); nrows=nrows, samples=samples)
end

println("\nsampled row access")
for (label, ix) in index_sets(indices)
    bench("raw[$label]", length(ix), () -> sampled_rows(raw, ix); nrows=nrows, samples=samples)
end
