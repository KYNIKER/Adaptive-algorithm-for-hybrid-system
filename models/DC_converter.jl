#using SparseArrays: sparsevec
#using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
#using ReachabilityAnalysis: add_dimension
using LazySets

function loadDC_Converter()
    d_aux = 0.0
    R_L = 0.3 # internal resistor
    R = 0.0 # load resistor

    L = 0.0 # inductance

    C_0 = 0.0

    x_1_min = 0.0
    x_1_max = 0.0
    x_2_min = 0.0
    x_2_max = 0.0


    x_1 = 1 # this is i_L(t)
    x_2 = 2 # this is v_0(t)

    dims = 2
    A = zeros(dims, dims)
    A_1 = [-R_L/L -1/L; 1/C_0 -1/(C_0*R)]
    A_2 = [0.0 1/L; -1/C_0 0.0]
    B = [1/L, 0.0]
    v_s = Zonotope(B, zeros(3, 1))

    S = Hyperrectangle([(x_1_max + x_1_min)/2, (x_2_max + x_2_min)/2], [(x_1_max - x_1_min)/2, (x_2_max - x_2_min)/2])
    locations = ModalLocation[]

    unsafe = [HPolyhedron([LazySets.HalfSpace([1.0, 0.0], x_1_min)]), HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], -x_1_max)]), HPolyhedron([LazySets.HalfSpace([0.0, 1.0], x_2_min)]), HPolyhedron([LazySets.HalfSpace([0.0, -1.0], -x_2_max)])]

    ModalHybridSystem(S, unsafe, locations), period
end

# The dc-dc boost converter of Karamanakos, Geyer and Manias, "Direct Voltage Control of DC–DC Boost Converters Using
# Enumeration-Based Model Predictive Control", IEEE Trans. Power Electron. 29(2), 2014, doi:10.1109/TPEL.2013.2256370,
# as a ModalHybridSystem (modal.jl). Returns the system and the decision period.
#
# The parameters and the safety box default to those of the DC-DC converter in Brorholt et al.'s "Shielded Learning for
# Hybrid Systems" code (Shared Code/DC-DC Converter.jl, taken from the UPPAAL model of the related "Proving and Improving
# Model Predictive Control of DC-DC Boost Converters"): a decision every 50 µs, i_L in [0, 4] A, v_o in [14.5, 15.5] V.
# There the load R varies in [30, 200] Ω; here it is a constant (default 73 Ω, their initial state), since a varying R
# makes the dynamics bilinear.
#
# Continuous state x = (i_L, v_o): inductor current and output voltage. Time is in units of time_unit seconds
# (microseconds by default, so the periods and step lengths stay well above the shields' absolute tolerances).
# The switch position u is chosen at the start of every period (the MPC's sampling interval T_s):
#   - continuous conduction, switch off (u = 0): dx/dt = A₁ x + B v_s, the diode conducts
#   - continuous conduction, switch on  (u = 1): dx/dt = (A₁ + A₂) x + B v_s
#   - discontinuous conduction (switch off, i_L = 0): dx/dt = A₃ x, the diode blocks and the current stays at 0
# with A₁ = [-R_L/L  -1/L; 1/C_o  -1/(C_o R)], A₂ = [0  1/L; -1/C_o  0], B = [1/L; 0], A₃ = [0  0; 0  -1/(C_o R)].
# Layers: switch off (BOOST_OFF) and on (BOOST_ON). Waiting keeps the switch where it is, the action turns it over.
# The off layer has the conduction modes as two locations: when the current reaches 0 with v_o >= v_s, an edge enters
# discontinuous conduction (the current is set to exactly 0), and once v_o drops to v_s the diode conducts again. The
# discontinuous location is only entered by that edge (starts = false), so every period starts in continuous conduction;
# with v_o >= v_s at i_L = 0 it is at the guard at once. The face i_L >= 0 of the state space is bordered by that guard.
#
# Safety: the current stays below i_max and the output voltage inside [v_min, v_max]. v_s_spread makes the input voltage
# an interval v_s ± v_s_spread.
const BOOST_OFF, BOOST_ON = 1, 2

function loadBoostConverterKaramanakos(; v_s=10.0, L=450e-6, R_L=0.3, C_o=220e-6, R=73.0, T_s=50e-6, time_unit=1e-6,
    i_max=4.0, v_min=14.5, v_max=15.5, v_s_spread=0.0)
    i_L, v_o = 1, 2
    n = 2
    e(i) = [j == i ? 1.0 : 0.0 for j in 1:n]
    keep_all = Matrix{Float64}(I, n, n)
    reset_current = diagm([0.0, 1.0])

    # The paper's matrices, in time units of time_unit
    A₁ = time_unit .* [-R_L/L -1/L; 1/C_o -1/(C_o*R)]
    A₂ = time_unit .* [0.0 1/L; -1/C_o 0.0]
    A₃ = time_unit .* [0.0 0.0; 0.0 -1/(C_o*R)]
    B = time_unit .* [1/L, 0.0]
    period = T_s / time_unit

    supply = Zonotope(B .* v_s, reshape(B .* v_s_spread, n, 1))
    no_input = Zonotope(zeros(n), zeros(n, 1))

    off_ccm, on, off_dcm = 1, 2, 3
    # The current reaches 0 while the diode would be reverse biased
    to_dcm = Edge(off_dcm, HPolyhedron([LazySets.HalfSpace(e(i_L), 0.0), LazySets.HalfSpace(-e(v_o), -v_s)]), reset_current, zeros(n))
    # The diode is forward biased again (with an interval v_s, as soon as it can be)
    to_ccm = Edge(off_ccm, HPolyhedron([LazySets.HalfSpace(e(v_o), v_s + v_s_spread)]), keep_all, zeros(n))
    turn_on = ModalAction(Action(1, nothing, keep_all, zeros(n)), on)
    turn_off = ModalAction(Action(2, nothing, keep_all, zeros(n)), off_ccm)

    locations = [
        ModalLocation(BOOST_OFF, nothing, A₁, supply, [to_dcm], [turn_on]),
        ModalLocation(BOOST_ON, nothing, A₁ + A₂, supply, Edge[], [turn_off]),
        ModalLocation(BOOST_OFF, nothing, A₃, no_input, [to_ccm], ModalAction[]; starts=false),
    ]

    unsafe = [HPolyhedron([LazySets.HalfSpace(-e(i_L), -i_max)]),
        HPolyhedron([LazySets.HalfSpace(e(v_o), v_min)]), HPolyhedron([LazySets.HalfSpace(-e(v_o), -v_max)])]

    # A margin around the safe box of more than a period's change, so the reach sets of safe cells stay on the grid
    v_high = v_s + v_s_spread
    di = (v_high + R_L * i_max) / L * T_s
    dv = (i_max + v_max / R) / C_o * T_s
    margin = 2 .* [di, dv]
    low_corner = [0.0, v_min - margin[2]]
    high_corner = [i_max + margin[1], v_max + margin[2]]
    S = Hyperrectangle((low_corner .+ high_corner) ./ 2, (high_corner .- low_corner) ./ 2)
    return ModalHybridSystem(S, unsafe, locations), period
end
