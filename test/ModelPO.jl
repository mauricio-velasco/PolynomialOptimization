"""
Checks that the constructor of ModelPOs throws errors when invoked with incore

"""
R=SOSRing(1,1,[3,3,3]) #Any SOS Ring would do for this test
modelPOObj = PolynomialOptimization.initializeModelPO(R,"SOS", "modelo1")
#The constraint type must be listed
@test_throws(ErrorException, PolynomialOptimization.initializeModelPO(R,"ERROR", "modelo1"))
#The solver must be allowed
@test_throws(ErrorException, PolynomialOptimization.initializeModelPO(R,"SOS", "modelo1", solverName = "ERROR"))
#There should be as many constraint types as gSOS
@test_throws(ErrorException, PolynomialOptimization.initializeModelPO(R,["SOS","SOS"], "modelo1"))
"""
Next we verify that the JuMP models produced by the ModelPO class are working.
Each solver should have its own set of tests

We begin with tests for the SCS solver
"""
R=SOSRing(1,1,[2,2,2]) #Any SOS Ring would do for this test
modelPOObj = PolynomialOptimization.initializeModelPO(R,"SOS", "modelo1", solverName = "SCS", solverAccuracy = 1e-10)
mObj = modelPOObj.modelObj
#verifies that a trivial model, optimizing 0 without constraints, gives optimal result
@test(PolynomialOptimization.solve(mObj)==:Optimal)

"""
Univariate Optimization of random polynomials

"""
#Simple instance of univariate optimization
#First set up a ring
R = SOSRing(1,1,4)
One = R.One
X=R.X
Y=R.Y
gSOS=R.gSOS
#Then set up a restriction polynomial
p1 = X(1)^8+34*X(1)^6+4*X(1)^2+5.4345
RP = p1-Y(1)-gSOS(1)
Obj = Y(1)
#Next we carry Out the optimization
modelPO1 = PolynomialOptimization.initializeModelPO(R,"SOS", "modelo1", solverAccuracy=1e-8)
addEqZeroConstraint(modelPO1, RP)
addObjective(modelPO1, "Max", Obj)
@test(PolynomialOptimization.solve(modelPO1.modelObj)==:Optimal)
@test_approx_eq_eps(5.4345, PolynomialOptimization.getobjectivevalue(modelPO1.modelObj), 1e-7)
