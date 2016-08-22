"""
Checks that univariate sum-of-squares polynomials are correctly implemented.
"""
d=rand(1:10)
R = SOSRing(1,0,d)
gSOS = R.gSOS
X = R.X
Z = R.Z
Mons = [X(1)^k for k in 0:d]
MatrixM = [Z(1,i,j) for i in 1:d+1, j in 1:d+1]
W = MatrixM*Mons
Res = ((W')*(Mons')')[1]
@test length((Res-gSOS(1)).terms) == 0

"""
Verifies that the initialized object has interoperable parts of the right types
"""
n=rand(1:5)
q=rand(1:3)
degreesVec=[rand(2:4), rand(2:4)]
# The most basic test should guarantee that the generic sums of squares are correct, explicitly
numericTypeOfRing = rand([Float32, Float64, Int32, Int64])
R = SOSRing(n,q,degreesVec, numericTypeOfRing)
One = R.One
X=R.X
Y=R.Y
Z=R.Z
gSOS=R.gSOS
d=38
for k in 1:n
  @test(typeof(X(k))==MPoly{numericTypeOfRing})
end
@test_throws(BoundsError , X(n+1))

for k in 1:q
  @test(typeof(Y(k))==MPoly{numericTypeOfRing})
end
@test_throws(BoundsError, Y(q+1))

x_1, y_1 = generators(MPoly{Float64}, :x_1, :y_1)
h= X(1)+Y(1) - MultiPoly.MPoly{Float64}(x_1 + y_1)
@test (length(h.terms) == 0)
#Verifies the symmetry and sizes of the matrix entries
for k in 1:length(degreesVec)
  matrixSize = binomial(n+degreesVec[k],n)
  for j in 1:matrixSize
      for i in 1:matrixSize
        @test(Z(k,i,j)==Z(k,j,i))
      end
  end
  @test_throws(BoundsError, Z(matrixSize+1,1,1))
end
@test_throws(BoundsError, gSOS(length(degreesVec)+1))
