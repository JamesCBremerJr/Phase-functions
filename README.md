This is a reference implementation of the algorithm described in 

      James Bremer "On the numerical solution of second order linear ordinary
      differential equations in the high-frequency regime," Applied
      and Computational Harmonic Analysis, 44 (2018), 312-349.

It allows for the solution of a large class of second order linear ordinary
differential equations in time independent of the frequency of oscillations
of their solutions.

Please note that this is an *unoptimized* implementation of the algorithm
meant princpally to demonstrate that the algorithm is independent of frequency.
There are a number of obvious inefficiencies.

Each of the files which comprise this package is briefly described below.  
More complete documentation can be found in the files themselves.    The principal routines 
for constructing nonoscillatory phase functions are located in the file kummer.f90.
Several example programs which make use of these routines are included in the package;
see laguerre_quad.f90, jacobi_quad.f90 and eval_legendre.f90.

------------------------------------------------------------------------------------------

1.  The file utils.f90 contains an assortment of utility routines for printing,
plotting and the like.

2.  The file chebyshev.f90 contains code for representing functions via piecewise
Chebyshev expansions.  

3.  The file tensor.f90 contains code for representing functions of two variables
as sums of tensor products of Chebyshev polynomials --- that is, in the form

   sum T_i(x) T_j(y)
    i,j

with T_n the Chebyshev polynomial of degree n.

4.  The file odesolve.f90 contains subroutines for solving quasilinear second order
ordinary differential equations --- i.e., those of the form

   y''(t) = f(t,y(t),y'(t))

---- via a fairly robust, but somewhat slow method.

5.  The file kummer.f90 contains code for constructing nonoscillatory phase
functions for second order linear ordinary differential equations of the form

   y''(t) + q(t) y(t) = 0     for all t in [a,b],

where q(t) is positive and smooth.

6.  The file laguerre_quad.f90 contains code for constructing Gauss-Laguerre
quadratures of large orders.

7.  The file jacobi_quad.f90 contains code for constructing Gauss-Jacobi quadrature
rules of large orders.

8.  The file eval_legendre.f90 contains code for evaluating Legendre functions
of the first and second kinds of large orders.  It depends on the module
kummer_legendre contained in the file kummer_legnedre.f90.
