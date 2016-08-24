# PolynomialOptimization

## A package for modeling and solving multivariate Sum-of-Squares optimization problems (a.k.a. SOS programs)

A Sum-of-squares optimization problem is one of the form

min_{y\in S}(c_1y_1 + ... + c_qy_q)

where S is the set of $(y_1,\dots, y_q) \in \mathbb{R}^q$ such that the polynomials g_j(x) given by g_j(x) = y_1 p_{j,1}(x)+ ... + y_q p_{j,q}(x) + p_{j,0}(x) are sums-of-squares of polynomials with real coefficients in the variables (x_1, ... , x_n)=:x for j=1,...,m.

In the above problem the constant c's and the polynomials p_{i,j}(x) for 1<=j<=m and 0<=i<=q are given and have real coefficients.

## Example:  
Use SOS programming to find the minimum value of the polynomial x^4 + x^3 - 2x^2 + 4.

To solve this problem we use the SOS program
max y_1  where y_1 is such that
g_1(x):= (x^4+x^3-2x^2+4) - y_1 is SOS
