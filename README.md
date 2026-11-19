# Adaptive-algorithm-for-hybrid-system

ReACT, an adaptive reachability analysis for hybrid systems, and shield synthesis built on it.

## Layout

```
src/ReACTShields.jl      the shield module: include("src/ReACTShields.jl"); using .ReACTShields
src/core/                shared by both parts: hybrid system types, discretization, the grid and cell covers
src/shields/             the four shields and their building blocks
models/                  shield models: bouncing ball, oil pump, cruise control
scripts/                 shield_sandbox.jl for experiments, setup.jl to install the packages
test/runtests.jl         regression tests for the shields
reachability/            the reachability benchmarks: ReACTed.jl, testsuite.jl, JuliaReach/, plotFuncs/ and their models
```

| Shield | File | For |
|---|---|---|
| `ReACTed_reachable_cell_t` | shields/pred_fast.jl | computes every cell (PredFast with `online=false`) |
| `ReACTed_reachable_cell_online_pred_fast` | shields/pred_fast.jl | only the cells that can become unsafe |
| `ReACTed_reachable_cell_adaptive_fast` | shields/adaptive_fast.jl | online, with cells refined where it helps |
| `ReACTed_reachable_cell_modal` | shields/modal.jl | systems with several locations and discrete modes |

Every shield returns `grid, allowed, lost, info`. `allowed` holds a bit per choice for each cell: bit 0 is waiting a
period, bit `i` the `i`-th action. A cell where nothing is allowed is dead. See `choice_bits` in src/shields/flat.jl.
In every shield, trying an action outside its guard does nothing, so there it is allowed exactly when waiting is. In a
cell partly inside the guard, the states outside it wait, so the action's successors include the waiting ones.
Each prints the time per phase and a short summary at the end; pass `verbose=false` to keep it quiet.

## Usage

Run everything from the repository root.

```
julia scripts/setup.jl                       # once, installs the packages
julia --project=. test/runtests.jl           # about a minute
julia --project=. reachability/testsuite.jl  # the reachability benchmarks, writes results/
```

`scripts/shield_sandbox.jl` loads the module and the models and runs a shield. It uses JLD2 to cache grids: set
`fresh_grid = true` the first time to compute the grid that is later loaded. With Revise, `includet` instead of
`include` picks up edits to src/ without restarting Julia.

The reachability side is not a module: `reachability/ReACTed.jl` and `reachability/ReACTedReachability.jl` define
helpers with the same names, and the benchmarks rely on including both in order.
