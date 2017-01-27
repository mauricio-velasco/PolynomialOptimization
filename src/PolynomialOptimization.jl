module PolynomialOptimization

using DataStructures
using MultiPoly
Base.promote_op{R}(::Base.MulFun, ::Type{MultiPoly.MPoly{R}}, ::Type{MultiPoly.MPoly{R}}) = MultiPoly.MPoly{R}
using JuMP
using SCS

export SOSRing, coefficientsInXV, initializeModelPO, solverStringFunc, addEqZeroConstraint, addObjective

include("SOSRing.jl")
include("ModelPO.jl")
include("ringsToModels.jl")

end
