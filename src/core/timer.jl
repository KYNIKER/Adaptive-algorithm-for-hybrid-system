mutable struct PhaseTimer
    start::UInt64
    last::UInt64
    names::Vector{String}
    seconds::Vector{Float64}
end

PhaseTimer() = (t = time_ns(); PhaseTimer(t, t, String[], Float64[]))

function lap!(t::PhaseTimer, name::String)
    now = time_ns()
    i = findfirst(==(name), t.names)
    if i === nothing
        push!(t.names, name)
        push!(t.seconds, 0.0)
        i = length(t.names)
    end
    t.seconds[i] += (now - t.last) / 1e9
    t.last = now
    return t
end

function report(t::PhaseTimer, lines::AbstractString...)
    phases = join(("$name $(round(s, digits=3)) s" for (name, s) in zip(t.names, t.seconds)), ", ")
    println("Total $(round((time_ns() - t.start) / 1e9, digits=3)) s: ", phases)
    foreach(println, lines)
end
