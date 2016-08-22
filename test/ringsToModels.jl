using Base.Test
using MultiPoly
using PolynomialOptimization

"""
Verifies that factoring works as expected on a
"""
n=2
q=2
d=3
R=SOSRing(n,q,d)
q = R.X(1)^2*Y(1)+R.X(1)^2*Z(1,1,1)
q2 = coefficientsInXV(R,q)
V = sum(map((k)->prod([R.X(j)^k[j] for j in 1:length(k)])*q2[k], keys(q2)))
@test length((V-q).terms) == 0

"""
Verifies that using the coefficientsInXV of a polynomial with non-linear terms throws an error
"""
q = R.X(1)^2*Y(1)^2
@test_throws(AssertionError, coefficientsInXV(R,q))
