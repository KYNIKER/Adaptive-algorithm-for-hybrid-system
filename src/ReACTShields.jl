"""
    ReACTShields

Shield synthesis for hybrid systems on a grid, built on the ReACT reachability analysis.

Shields (each returns `grid, allowed, lost, info`, see `choice_bits` in shields/flat.jl):
  - `ReACTed_reachable_cell_t`                : computes every cell (shields/pred_fast.jl, with online = false)
  - `ReACTed_reachable_cell_online_pred_fast` : online, only the cells that can become unsafe (shields/pred_fast.jl)
  - `ReACTed_reachable_cell_adaptive_fast`    : online and multi-resolution (shields/adaptive_fast.jl)
  - `ReACTed_reachable_cell_modal`            : several locations and discrete modes (shields/modal.jl)

Use it from a script with

    include("src/ReACTShields.jl"); using .ReACTShields

or, to have edits picked up without restarting Julia, `using Revise; includet("src/ReACTShields.jl"); using .ReACTShields`.
"""
module ReACTShields

using LazySets, LinearAlgebra, ReachabilityAnalysis, Dates, Plots
using SparseArrays: sparsevec

# Shared with the reachability benchmarks in reachability/
include("core/discretize.jl")
include("core/utilities.jl")
# The grid, cells and covers of zonotopes by cells
include("core/partition.jl")
include("core/timer.jl")          # PhaseTimer, the time per phase each shield reports

# Order matters: later files use the types and helpers of earlier ones.
include("shields/guards.jl")        # ReACT_guards_t and the set predicates (allunder, someinside, ...)
include("shields/propagate.jl")     # reach set of a cell over one period, with jumps
include("shields/flat.jl")          # successor store, the shield state and its rules, choice_bits
include("shields/pred_fast.jl")
include("shields/adaptive_fast.jl")
include("shields/modal.jl")

include("../models/bouncingBall.jl")
include("../models/oilpump.jl")
include("../models/cruise.jl")

export ReACTed_reachable_cell_t, ReACTed_reachable_cell_online_pred_fast, ReACTed_reachable_cell_adaptive_fast, ReACTed_reachable_cell_modal
export Grid, Edge, Action, EuclideanHybridSystem, ModalHybridSystem, ModalLocation, ModalAction
export loadBouncingBallShielded, loadBouncingBallShieldedWithInput, loadBouncingBallNoEdges
export loadOilPumpModal, plot_oil_pump, OILPUMP_OFF_FREE, OILPUMP_OFF_LOCKED, OILPUMP_ON_FREE, OILPUMP_ON_LOCKED, OILPUMP_PIECES
export loadCruiseControl, loadCruiseControlModal, loadCruiseControlModal_t, CRUISE_DECELERATE, CRUISE_NEUTRAL, CRUISE_ACCELERATE

end
