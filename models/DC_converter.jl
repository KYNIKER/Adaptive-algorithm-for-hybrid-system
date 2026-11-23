#using SparseArrays: sparsevec
#using ReachabilityAnalysis.ReachabilityBase.Arrays: SingleEntryVector
#using ReachabilityAnalysis: add_dimension
using LazySets

include("../src/core/utilities.jl")

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
    A_1 = [-R_L/L, -1/L; 1/C_0, -1/(C_0*R)]
    A_2 = [0.0, 1/L; -1/C_0, 0.0]
    B = [1/L, 0.0]
    v_s = Zonotope(B, zeros(3, 1))

    S = Hyperrectangle([(x_1_max + x_1_min)/2, (x_2_max + x_2_min)/2], [(x_1_max - x_1_min)/2, (x_2_max - x_2_min)/2])
    locations = ModalLocation[]

    unsafe = [HPolyhedron([LazySets.HalfSpace([1.0, 0.0], x_1_min)]), HPolyhedron([LazySets.HalfSpace([-1.0, 0.0], -x_1_max)]), HPolyhedron([LazySets.HalfSpace([0.0, 1.0], x_2_min)]), HPolyhedron([LazySets.HalfSpace([0.0, -1.0], -x_2_max)])]

    ModalHybridSystem(S, unsafe, locations), period
end