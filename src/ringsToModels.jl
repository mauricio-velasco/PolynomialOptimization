"""
coefficientsInXV(multivarPoly,R::SOSRing)

Given multivarPoly in k[X,Y,Z] returns its vector of coefficients in k [Y,Z] [X],
if these are affine LINEAR in the Y's and Z's, else returns ERROR.
"""
function coefficientsInXV(R::SOSRing, multivarPoly)
  Res=Dict()
  n=R.n
  for (k,v) in multivarPoly.terms
      a = k[1:n]
      b = k[n+1:end]
      if sum(b)>=2
        throw(AssertionError("Some coefficient of a monomial in Xs is NOT an affine linear function of the remaining variables."))
      end
      haskey(Res,a) ? Res[a][b] = v : Res[a] = Dict(b=>v)
  end
  Res2 = Dict()
  for mon in keys(Res)
    Res2[mon] = MultiPoly.MPoly{Float64}(OrderedDict(Res[mon]),R.LateSymbols)
  end
  return Res2
end

#=
function addLinearConstraint(RP, R::SOSRing)
  #Adds the constraints encoded by a polynomial RP which is linear in Y and the Z_i
  Res = getLinearConstraints(RP,R)
  stringRes = []
  for mon in keys(Res)
    st = string(Res[mon])*" == 0"
    push!(stringRes, st)
  end
  for st in stringRes
    eval(parse("@constraint(m,"*st*")"))
  end
end

function addObjective(Obj, MaxOrMin="Max")
  if MaxOrMin == "Min"
    eval(parse("@objective(m, Min,"*string(Obj)*")"))
  elseif MaxOrMin == "Max"
    eval(parse("@objective(m, Max,"*string(Obj)*")"))
  else
    print("Please specify Max or Min")
  end
end

function addTraceRest(R::SOSRing)
  #If we are only interested in a spectrahedral cone we need to choose a grading.
  #The most natural is to normalize de trace. The most common case here is that the objective be the 0 function.
  numDegrees = length(R.degreesVec)
  degreesVec = R.degreesVec
  Res = []
  for k in 1:numDegrees
    for i in 1:degreesVec[k]
      push!(Res, R.Z(k,i,i))
    end
  end
  RP= sum(Res)
  eval(parse("@constraint(m,"*string(RP)*"==1.0)"))
  eval(parse("@objective(m, Max,0.0)"))
end

function evaluateMonomAtOpt(exponentVector, R::SOSRing)
  #Part of a function
  n = R.n
  a = exponentVector[1:n]
  b = exponentVector[n+1:end]
  Coeff = 1
  for j in 1:length(b)
    VarValue = eval(parse("getvalue("*string(R.VarSymbols[n+j])*")"))
    Coeff = Coeff * VarValue^b[j]
  end
  XSymbols = R.VarSymbols[1:n]
  return MPoly{Float64}(OrderedDict(a => Coeff),XSymbols)
end

function evaluateAtOpt(multivarPoly,R::SOSRing)
  #Evaluates a polynomial at the optimum. Generally the result is a polynomial involving only the X variables.
  #TODO: Needs to check the status of the solver before proceeding.
  Res=[]
  Trms = terms(multivarPoly)
  for (k,v) in Trms
      push!(Res, v * evaluateMonomAtOpt(k,R))
  end
  return sum(Res)
end

function biggestCoeff(multivarPoly)
  Res = 0
  for (k,v) in multivarPoly.terms
    if abs(v)>Res
      Res=abs(v)
    end
  end
  return Res
end

function relevantTerms(multivarPoly, R::SOSRing, Tol=0.0001)
  #TODO: Write a function which returns the summands in the sum-of-squares
  #TODO: Write function which rounds a polynomial, producing one with fewer monomials, if possible
end
=#
