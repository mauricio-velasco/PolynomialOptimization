module PolynomialOptimization

using DataStructures
using MultiPoly
using JuMP
using SCS

export SOSRing, coefficientsInXV, initializeModelPO, solverStringFunc, addEqZeroConstraint, addObjective

include("SOSRing.jl")
include("ModelPO.jl")
include("ringsToModels.jl")

end
