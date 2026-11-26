"""
    ReACTShields

Shield synthesis for hybrid systems on a grid, built on the ReACT reachability analysis.

Shields (each returns `grid, allowed, lost, info`, see `choice_bits` in shields/flat.jl):
  - `ReACTed_reachable_cell_modal`            : several locations and discrete modes (shields/modal.jl), offline or online
  - `ReACTed_reachable_cell_modal_adaptive`   : the same, multi-resolution (shields/modal_adaptive.jl)
  - `ReACTed_reachable_cell_modal_simple`     : the same method without optimizations or clock axes (shields/modal_simple.jl)
A system without discrete modes (an `EuclideanHybridSystem`, such as the bouncing ball) runs through them as
`ModalHybridSystem(system)`.

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
include("core/partition.jl")
include("core/timer.jl")

# Order matters: later files use the types and helpers of earlier ones.
include("shields/guards.jl")
include("shields/propagate.jl")
include("shields/flat.jl")
include("shields/modal_shared.jl")
include("shields/modal.jl")
include("shields/modal_adaptive.jl")
include("shields/modal_simple.jl")

include("../models/bouncingBall.jl")
include("../models/oilpump.jl")
include("../models/cruise.jl")
include("../models/DC_converter.jl")

export ReACTed_reachable_cell_modal, ReACTed_reachable_cell_modal_adaptive, ReACTed_reachable_cell_modal_simple
export Grid, Edge, Action, EuclideanHybridSystem, ModalHybridSystem, ModalLocation, ModalAction
export loadBouncingBallShielded, loadBouncingBallShieldedWithInput, loadBouncingBallNoEdges
export loadOilPumpModal, plot_oil_pump, OILPUMP_OFF_FREE, OILPUMP_OFF_LOCKED, OILPUMP_ON_FREE, OILPUMP_ON_LOCKED, OILPUMP_PIECES
export loadCruiseControl, loadCruiseControlModal, loadCruiseControlModal_t, plot_cruise, CRUISE_DECELERATE, CRUISE_NEUTRAL, CRUISE_ACCELERATE
export loadBoostConverterKaramanakos, BOOST_OFF, BOOST_ON

end
