# Times a shield on a small model and profiles where its time goes, in one run. Run from the repository root:
#
#     julia --project=. scripts/profile_shields.jl <model> [online=true] [allocs=false]
#
# Models (small on purpose, so a profile takes a few minutes at most):
#   cruise_t       loadCruiseControlModal_t, sensor range 50, g = 1, δ 0.25-1, act_needs_safe_flow = false (modal)
#   cruise         loadCruiseControlModal, sensor range 50, g = 1, δ 0.25-1 (modal)
#   oilpump        loadOilPumpModal, period 0.2, g = 0.2, δ 0.2-0.4 (modal)
#   ball_pred      bouncing ball with input, g = 0.03, δ 0.001-0.032 (pred_fast)
#   ball_adaptive  bouncing ball with input, g0 = 0.06, L = 1, δ 0.001-0.032 (adaptive_fast)
# The ball's adaptive shield is always online.
#
# It compiles the code with a warm-up on a coarser grid (twice the granularity), then does one run under @timed and the
# CPU profiler together: @timed counts the bytes and the GC time exactly, the profiler samples where the time goes, and
# its overhead (a sample every 2 ms) adds a few percent to the time. It prints the time, bytes, GC time and the GC share
# of the samples, and writes results/profile/<model>_cpu.txt: the CPU profile, flat, by sample count (inclusive count and
# self "Overhead" per line), with C frames, so the garbage collector shows up under the code that triggered it.
# With allocs=true it then does a second run under the allocation profiler and writes results/profile/<model>_alloc.txt:
# a 1% sample of the allocations by the innermost frame in src/ or models/, with their types. That profiler slows the
# code that allocates much more than the rest, so it gets a run of its own; it is for finding where allocations come
# from, not for timing.
# Baselines (warm, separate runs) of 2026-11-17, before the optimizations of 2026-11-18: cruise_t offline 29 s and
# 16.1 GiB, oil pump 3.6 s and 2.7 GiB, ball_pred 2.7 s and 2.5 GiB, ball_adaptive 4.1 s and 3.5 GiB.
using LazySets, ReachabilityAnalysis, Profile, Printf
const REPO = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "src", "ReACTShields.jl"))
using .ReACTShields

const MODELS = ("cruise_t", "cruise", "oilpump", "ball_pred", "ball_adaptive")
isempty(ARGS) && error("give a model: $(join(MODELS, ", "))")
model = ARGS[1]
model in MODELS || error("unknown model $model; the models are $(join(MODELS, ", "))")
opts = Dict(split(a, "="; limit=2) for a in ARGS[2:end])
online = parse(Bool, get(opts, "online", "true"))
allocs = parse(Bool, get(opts, "allocs", "false"))

# scale = 2 is the coarse warm-up
if startswith(model, "ball")
    sys, T = loadBouncingBallShieldedWithInput()
    synthesize(scale=1) = model == "ball_pred" ? ReACTed_reachable_cell_online_pred_fast(sys, T, 0.03scale, 0.001, 0.032; online, verbose=false) :
                          ReACTed_reachable_cell_adaptive_fast(sys, T, 0.06scale, 1, 0.001, 0.032; patience=0, verbose=false)
else
    sys, T, g, δ⁻, δ⁺, ansf = model == "cruise_t" ? (loadCruiseControlModal_t(; sensor_range=50.0)..., 1.0, 0.25, 1.0, false) :
                             model == "cruise" ? (redirect_stdout(() -> loadCruiseControlModal(; sensor_range=50.0), devnull)..., 1.0, 0.25, 1.0, true) :
                             (loadOilPumpModal(; period=0.2)..., 0.2, 0.2, 0.4, true)
    synthesize(scale=1) = ReACTed_reachable_cell_modal(sys, T, scale * g, δ⁻, δ⁺; online, act_needs_safe_flow=ansf, verbose=false)
end

synthesize(2) # warm-up, so the run leaves out compilation
GC.gc()
outdir = joinpath(REPO, "results", "profile")
mkpath(outdir)
prefix = joinpath(outdir, model)

Profile.clear()
Profile.init(n=10^7, delay=0.002)
st = @timed @profile synthesize()
allowed = st.value[2]
open("$(prefix)_cpu.txt", "w") do io
    Profile.print(IOContext(io, :displaysize => (100000, 400)); format=:flat, sortedby=:count, mincount=20, C=true)
end
# Share of the samples in the garbage collector, from the backtraces
data = Profile.fetch(include_meta=false) # backtraces separated by a single 0
lidict = Profile.getdict(data)
gc_frame(ip) = any(fr -> occursin("gc_collect", string(fr.func)), lidict[ip])
let in_gc = false, n = 0, n_gc = 0
    for ip in data
        if ip == 0
            n += 1; n_gc += in_gc; in_gc = false
        elseif !in_gc && haskey(lidict, ip)
            in_gc = gc_frame(ip)
        end
    end
    global gc_share = n == 0 ? 0.0 : n_gc / n
end
@printf("%s online=%s: %.2f s (with the profiler), %.2f GiB allocated, GC %.2f s (%.0f%%), %.0f%% of samples in GC, dead %d of %d
",
    model, online, st.time, st.bytes / 2^30, st.gctime, 100 * st.gctime / st.time, 100 * gc_share, count(==(0), allowed), length(allowed))
println("wrote $(prefix)_cpu.txt")
allocs || exit()

# Allocation profile: a 1% sample, attributed to the innermost frame in src/ or models/
in_repo(file) = (f = lowercase(normpath(string(file))); any(d -> startswith(f, lowercase(joinpath(REPO, d))), ("src", "models")))
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 0.01 synthesize()
res = Profile.Allocs.fetch()
bytes = Dict{String,Int}(); counts = Dict{String,Int}(); types = Dict{String,Dict{String,Int}}()
for a in res.allocs
    fr = findfirst(f -> in_repo(f.file), a.stacktrace)
    key = fr === nothing ? "outside src/ and models/" : (f = a.stacktrace[fr]; string(f.func, " @ ", basename(string(f.file)), ":", f.line))
    bytes[key] = get(bytes, key, 0) + a.size
    counts[key] = get(counts, key, 0) + 1
    t = get!(types, key, Dict{String,Int}())
    t[string(a.type)] = get(t, string(a.type), 0) + 1
end
total_b = sum(values(bytes)); total_c = sum(values(counts))
open("$(prefix)_alloc.txt", "w") do io
    println(io, "sampled allocations: $total_c, bytes $(round(total_b / 2^20, digits=1)) MiB (1% sample)")
    for (k, v) in sort(collect(counts), by=last, rev=true)[1:min(40, end)]
        tt = join(("$(first(x))×$(last(x))" for x in sort(collect(types[k]), by=last, rev=true)[1:min(3, end)]), ", ")
        @printf(io, "%5.1f%% of allocs, %5.1f%% of bytes  %s  [%s]\n", 100v / total_c, 100bytes[k] / total_b, k, first(tt, 160))
    end
end
println("wrote $(prefix)_alloc.txt")
