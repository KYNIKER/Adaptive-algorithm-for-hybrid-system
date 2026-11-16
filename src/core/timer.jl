# Wall-clock time of the phases of a shield synthesis, for the summary each shield prints at the end.
# Each lap! adds the time since the previous lap (or since the timer was made) to a named phase, so a loop that does
# reachability then shield work just calls lap!(t, "Reachability") and lap!(t, "Shield") in turn.
mutable struct PhaseTimer
    start::UInt64
    last::UInt64
    names::Vector{String} # in the order the phases first appeared
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

# Prints the total time since the timer was made, the time of each phase and then each line of `lines`.
function report(t::PhaseTimer, lines::AbstractString...)
    phases = join(("$name $(round(s, digits=3)) s" for (name, s) in zip(t.names, t.seconds)), ", ")
    println("Total $(round((time_ns() - t.start) / 1e9, digits=3)) s: ", phases)
    foreach(println, lines)
end
