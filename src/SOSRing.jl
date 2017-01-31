"""
SOSRing

This is the type used for writing polynomial optimization
problems with sums of squares. Any instance of SOSRing contains

n => number of variables for the polynomials
q => number of optimization variables
degreesVec => vector of half-degrees of the sums of squares in the n variables above
One => function returning the unit of the ring
X => function returning variables X(1) ... X(n) for building polynomials
Y => function returning optimization variables Y(1) ... Y(q)
Z => function returning the matrix variables Z(k,i,j) (k-th symmetric matrix row i column j)
gSOS => function returning the sums of squares g(1) ... g(length(degreesVec))
VarSymbols => Array with ALL variables on our ring
LateSymbols => Array with ALL symbols except the X variables

Examples:

julia> SOSRingObj(2,1,[2,3]) #Ring with variables x_1, x_2, y_1 and sums of squares of degree 4 and 6 in the x's

julia> SOSRingObj(2,1,3) #Ring with variables x_1, x_2, y_1 and a sum of squares of degree 6 in the x's

"""
type SOSRing
  n
  q
  degreesVec
  One
  X
  Y
  Z
  gSOS
  VarSymbols
  LateSymbols
  modelName::AbstractString
end

function SOSRing(n,q,degreesVec::Array{Int}, numericTypeOfRing = Float64, modelName = "m")
  #Most often we know only how many variables n, constraints variables q and SOS variables
  #are needed but the constructor requires other input
  One, X, Y, Z, G, VarSymbols, LateSymbols = realSOSRing(n,q,degreesVec, numericTypeOfRing)
  return SOSRing(n,q,degreesVec,One, X,Y,Z,G,VarSymbols, LateSymbols, modelName)
end

function SOSRing(n,q,d::Int, numericTypeOfRing = Float64, modelName = "m")
  degreesVec=[d]
  One, X, Y, Z, G, VarSymbols, LateSymbols = realSOSRing(n,q,degreesVec, numericTypeOfRing)
  return SOSRing(n,q,degreesVec,One, X,Y,Z,G,VarSymbols,LateSymbols, modelName)
end

"""
exponentVctsAddOneMoreVar(monomialList, atMostd)

Given the array monomialList of exponent vectors of length k with sum <= atMostd
returns the array of all exponent vectors of length k+1 with sum <= atMostd
whose first k components were in some element of monomialList
"""
function exponentVctsAddOneMoreVar(monomialList, atMostd)
  Result = []
  for mon in monomialList
    t=sum(mon)
    for k in [0:atMostd-t;]
      monNew = copy(mon)
      push!(monNew,k)
      push!(Result,monNew)
    end
  end
  return Result
end

"""
monomialExponents(n,d)

Computes all exponent vectors of monomials of degree at most d in n variables
"""
function monomialExponents(n,d)
  Monoms=Array[BigInt[]]
  for k in [1:n;]
    Monoms=exponentVctsAddOneMoreVar(Monoms,d)
  end
  return Monoms
end

"""
intArray2MonomFn(Index2VariableFunc, indices)

Returns a monomial from the vector of exponent indices. The corresponding monomial lives in a
ring with variables given by the function Index2VariableFunc
"""
function intArray2MonomFn(Index2VariableFunc, indices)
  Res = Index2VariableFunc(1)^(0)
  for k in [1:length(indices);]
      Res = Res * (Index2VariableFunc(k))^(indices[k])
    end
  return Res
end

"""
ZVariables(n,degreesVec)

Returns two arrays of size length(degreesVec):
One containing the variables of all symmetric matrices encoding the gsos of degrees 2*d: d in degreesVec
and another array containing their exponent vectors (which will be appended to that of X's and Y's to
obtain the total degree in the final ring).
"""
function ZVariables(n,degreesVec)
  Zs = Symbol[]
  resIndices =[]
  for k = 1:length(degreesVec)
    d = degreesVec[k]
    binomSize = binomial(d+n,n)
    Zindices = Array[]
    dimSim2 = BigInt(binomSize*(binomSize+1)/2)
    counterEach = 1
    for i in [1:binomSize;]
      for j in [1:binomSize;]
        if i<=j
          baseV = vec(zeros(BigInt, 1, dimSim2))#This is the number of variables in a symmetric binomSize x binomSize matrix
          baseV[counterEach] = BigInt(1)
          counterEach = counterEach +1
          push!(Zindices, baseV)
          push!(Zs, Symbol("Z$k[$i,$j]"))
        end
      end
    end
    push!(resIndices,Zindices)
  end
  return Zs,resIndices
end

"""
startingIndex(n,q,k,degreesVec)

Returns the index of the variable Z(k,1,1) in the exponent vector
of all variables of our ring.
"""
function startingIndex(n,q,k,degreesVec)
  index = n+q
  for l in 1:k-1
    binomSize = binomial(n+degreesVec[l],degreesVec[l])
    index += BigInt(binomSize*(binomSize+1)/2)
  end
  return index
end

"""
realSOSRing(n, q, degreesVec, numericTypeOfRing = Float64)

Does the main work of the constructor of an SOSRing.
Returns functions XV,YV,g which access the variables X(1)...X(n), Y(1)...Y(q) and the sums
of squares g(i)=mon^tZ_i mon and the symbols corresponding to Variables X's, Y's and Z's.
As well as LateSymbols (which are the non-X variables).
"""
function realSOSRing(n, q, degreesVec, numericTypeOfRing = Float64)
  #TODO: Possible improvement. The function gSOS recomputes polynomials every time it is called. muss es sein?
  Xs = [Symbol("x_$i") for i in 1:n]
  Ys = [Symbol("y_$i") for i in 1:q]
  ZsArray, AllZIndices = ZVariables(n,degreesVec)
  Zs = vcat(ZsArray)
  VarSymbols  = vcat(Xs,Ys,Zs)#This array contains the symbolic display of variables
  LateSymbols = vcat(Ys,Zs)
  #totSize will become the total number of variables x_i,y_j,Z[j,k]
  totSize = n+q
  for d in degreesVec
    binomSize = binomial(d+n,d)
    dimSim2 = BigInt(binomSize*(binomSize+1)/2)
    totSize = totSize + dimSim2
  end
  T = numericTypeOfRing
  #Unit of the ring function
  function One()
    baseV = vec(zeros(BigInt, 1, totSize))
    return MPoly{T}(OrderedDict(baseV=> 1.0),VarSymbols)
  end
  #function which calls the i-th X variable (todo: Overwrite with metaprogramming)
  function XV(i)
    if i<=n && 0<i
      baseV = vec(zeros(BigInt, 1, totSize))
      baseV[i] = 1
      return MPoly{T}(OrderedDict(baseV=> 1.0),VarSymbols)
    else
      throw(BoundsError(["Index out of range"]))
    end
  end
  #function which calls the i-th y variable
  function YV(i)
    if i<=q && 0<i
      baseV = vec(zeros(BigInt, 1, totSize))
      baseV[i+n] = 1
      return MPoly{T}(OrderedDict(baseV=> 1.0),VarSymbols)
    else
      throw(BoundsError(["Index out of range."]))
    end
  end
  #Functions calling variables (i,j) in the k-th symmetric matrix
  function ZV(k,i,j)
    stIndex = startingIndex(n,q,k,degreesVec)
    binomSize = binomial(degreesVec[k]+n,degreesVec[k])
    if i<=binomSize && 0<i && 0<j && j<=binomSize
      baseV = vec(zeros(Int64, 1, totSize))
      b=binomSize
      if i<=j
        t= Int64(((b+1)*b)/2 - (b+2-i)*(b+1-i)/2+(j-i+1))
      else
        t= Int64(((b+1)*b)/2 - (b+2-j)*(b+1-j)/2+(i-j+1))
      end
      baseV[t+stIndex] = 1
      return MPoly{T}(OrderedDict(baseV=> 1.0),VarSymbols)
    else
      throw(BoundsError(["Index out of range."]))
    end
  end
  #function which builds the i-th generic sum of squares
  function gV(k)
    d=degreesVec[k]
    binomSize = binomial(d+n,d)
    ZMatrix = MultiPoly.MPoly{T}[ZV(k,i,j) for i in 1:binomSize, j in 1:binomSize]
    Mons = monomialExponents(n,d)
    MonFns = MultiPoly.MPoly{T}[]
    for index in Mons
      push!(MonFns, intArray2MonomFn(XV,index))
    end
    Res = MonFns'*ZMatrix
    SOSPolynomial = Res*MonFns
    GenericSOSP=SOSPolynomial[1]
    return GenericSOSP
  end

return One,XV,YV,ZV,gV, VarSymbols, LateSymbols
end
