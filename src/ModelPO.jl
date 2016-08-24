"""
ALLOWED GLOBAL PARAMETERS
"""
global allowedConstraints = ["SOS"]
# "SOS" means that the undelying symmetric matrix will be required to be PSD.
global allowedSolvers = ["SCS"]
"""
ModelPO

Each polynomial optimization problem is built in some ring R::SOSRing
Each of the gSOS of R determines a constraint in constraintTypes::Array{AbstractString}
which must be an element of the global variable allowedConstraints.

The corresponding JuMP model has name modelName::AbstractString
and is accessible directly via modelObj. The solver info is encoded in the solverString field.

Finally, for printing, the steps of the construction are recorded, as strings in
modelDescription::Array{AbstractString}
"""
type ModelPO
  R::SOSRing
  constraintTypes::Array{AbstractString}
  modelName::AbstractString
  modelObj
  modelDescription::Array{AbstractString}
  solverString
end
"""
solverString

Given a solverName and a numerical solverAccuracy returns
the description of the solver to be used by JuMP
"""
function solverStringFunc(solverName::AbstractString, solverAccuracy)
  if !(solverName in allowedSolvers) throw(ErrorException("Solver not recognized")) end
  stringAccuracy = string(solverAccuracy)
  if solverName == "SCS"
    solverString = "solver = SCSSolver(verbose = 0, eps = $stringAccuracy)"
  end
  return solverString
end
"""
initializeModelPO(R::SOSRing, constraintTypes::Array{AbstractString}, modelName::AbstractString, modelDescription = ["SOS"], solverName ="SCS", solverAccuracy = 1e-6)

This is the main function for constructing ModelPO objects.
"""
function initializeModelPO(R::SOSRing, constraintTypes::Array{ASCIIString}, modelName::AbstractString ; modelDescription = ["SOS"], solverName ="SCS", solverAccuracy = 1e-6)
  #First, we define the JuMP model
  solverString = solverStringFunc(solverName, solverAccuracy)
  eval(parse("global $modelName = Model($solverString)"))
  #We add the optimization variables
  for k in 1:R.q
    eval(parse("@variable($modelName,y_$k)"))
  end
  #Next we add the constraints, as specified in constraintTypes
  #Verify that the constraintTypes coincides with the constraints
  if length(constraintTypes) != length(R.degreesVec)
    throw(ErrorException("The constraintTypes vector an the gSOSs of R have a different number of components"))
  end
  for k in 1:length(constraintTypes)
    #Usual sums of squares, we require the corresponding symmetric matrix to be PSD
    if constraintTypes[k] == "SOS"
      d = R.degreesVec[k]
      binomSize = binomial(R.n + d,d)
      eval(parse("@variable($modelName,Z$k[1:$binomSize,1:$binomSize],SDP)"))
    end
  end
  eval(parse("modelVariable = PolynomialOptimization.$modelName"))
  return ModelPO(R, constraintTypes, modelName, modelVariable, modelDescription, solverString)
end

"""
Simpler constructor for ModelPO objects

Examples:
modelPOObj = initializeModelPO(R, "SOS", "modelo1")

"""
function initializeModelPO(R::SOSRing, constraintTypes::AbstractString, modelName::AbstractString ; modelDescription = ["SOS"], solverName ="SCS", solverAccuracy = 1e-6)
  if !(constraintTypes in allowedConstraints) throw(ErrorException("Constraint type not recognized")) end
  if constraintTypes == "SOS"
    constraintTypesVector = ASCIIString[constraintTypes for k in 1:length(R.degreesVec)]
    return initializeModelPO(R, constraintTypesVector, modelName, modelDescription = modelDescription,solverName = solverName, solverAccuracy = solverAccuracy)
  else
    throw(ErrorException("Unrecognized constraint type"))
  end
end

"""
"""
function addEqZeroConstraint(mod::ModelPO, RP::MultiPoly.MPoly)
  #Adds the constraints encoded by a polynomial RP which is linear in Y and the Z_i
  Res =  coefficientsInXV(mod.R, RP)
  modelName = mod.modelName
  stringRes = []
  for mon in keys(Res)
    st = string(Res[mon])*" == 0"
    push!(stringRes, st)
  end
  for st in stringRes
    eval(parse("@constraint($modelName , $st )"))
  end
end

function addObjective(mod::ModelPO, maxormin, Obj)
  if !(maxormin in ["Min", "Max"]) throw(ErrorException("Only Max or Min allowed on objective function")) end
  objectiveStr = string(Obj)
  modelName = mod.modelName
  eval(parse("@objective($modelName, $maxormin, $objectiveStr)"))
end
