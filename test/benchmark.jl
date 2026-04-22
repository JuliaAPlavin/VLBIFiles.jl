# julia --project=test test/benchmark.jl
using VLBIFiles

p = "/Users/aplavin/work/galactic_scatter/2005+403/data/archival/vlba/UG002/VLBA_UG002O_ug002o_BIN0_SRC0_0_180821T142741.idifits"
# p = joinpath(@__DIR__, "data/BL146_1.fits")
uvf = VLBI.load(VLBI.UVData, p)

t = @elapsed raw = VLBI.read_data_raw(uvf)
nrows = length(raw)
println("read_data_raw: $(round(t, digits=2))s  ($nrows rows)")

function report(label, n, t)
    us = t / n * 1e6
    s = t / n * nrows
    println("  $(rpad(label, 28)) $(lpad(round(us, digits=2), 10)) µs/elem  $(lpad(round(s, digits=2), 8)) s/file")
end

n = 10000

report("collect SOURCE (Int32)", nrows, @elapsed collect(raw.SOURCE))
report("collect DATE (Float64)", nrows, @elapsed collect(raw.DATE))

report("SOURCE[rand] ×$n", n, @elapsed begin for i in rand(1:nrows, n); raw.SOURCE[i]; end end)
report("SOURCE[sorted] ×$n", n, @elapsed begin for i in sort(rand(1:nrows, n)); raw.SOURCE[i]; end end)
report("SOURCE[rand×$n] ", n, @elapsed raw.SOURCE[rand(1:nrows, n)])
report("SOURCE[sorted×$n]", n, @elapsed raw.SOURCE[sort(rand(1:nrows, n))])

report("FLUX[rand] ×$n", n, @elapsed begin for i in rand(1:nrows, n); raw.FLUX[i]; end end)
report("FLUX[sorted] ×$n", n, @elapsed begin for i in sort(rand(1:nrows, n)); raw.FLUX[i]; end end)
report("FLUX[rand×$n] ", n, @elapsed raw.FLUX[rand(1:nrows, n)])
report("FLUX[sorted×$n]", n, @elapsed raw.FLUX[sort(rand(1:nrows, n))])

i0 = rand(1:max(nrows-n, 1))
report("row iteration ×$n", min(n, nrows), @elapsed begin for i in i0:min(i0+n-1, nrows); raw[i]; end end)
