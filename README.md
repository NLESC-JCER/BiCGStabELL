# BiCGStabELL
Krylov subspace iterative solvers for linear systems
bicgstabell() implements an iterative solver for non-symmetric linear
   operators, using the algorithm described in:

      Gerard L. G. Sleijpen and Diederik R. Fokkema, "BiCGSTAB(L) for
      linear equations involving unsymmetric matrices with complex
      spectrum," Electronic Trans. on Numerical Analysis 1, 11-32
      (1993).

   and also:

      Gerard L.G. Sleijpen, Henk A. van der Vorst, and Diederik
      R. Fokkema, " BiCGstab(L) and other Hybrid Bi-CG Methods,"
      Numerical Algorithms 7, 75-109 (1994).

   This is a generalization of the stabilized biconjugate-gradient
   (BiCGSTAB) algorithm proposed by van der Vorst (and described
   in the book _Templates for the Solution of Linear Systems_ by
   Barrett et al.)  BiCGSTAB(1) is equivalent to BiCGSTAB, and
   BiCGSTAB(2) is a slightly more efficient version of the BiCGSTAB2
   algorithm by Gutknecht, while BiCGSTAB(L>2) is a further
   generalization.

   The reason that we use this generalization of BiCGSTAB is that the
   BiCGSTAB(1) algorithm was observed by Sleijpen and Fokkema to have
   poor (or even failing) convergence when the linear operator has
   near-pure imaginary eigenvalues.  This is precisely the case for
   our problem (the eigenvalues of the timestep operator are i*omega),
   and we observed precisely such stagnation of convergence.  The
   BiCGSTAB(2) algorithm was reported to fix most such convergence
   problems, and indeed L > 1 seems to converge well for us. 

Other variations to explore:

   G. L. G. Sleijpen and H. A. van der Vorst, "Reliable updated
   residuals in hybrid Bi-CG methods," Computing 56 (2), 141-163
   (1996).

   G. L. G. Sleijpen and H. A. van der Vorst, "Maintaining convergence
   properties of BiCGstab methods in finite precision arithmetic,"
   Numerical Algorithms 10, 203-223 (1995).

   See also code on Sleijpen's web page:
                 http://www.math.uu.nl/people/sleijpen/
                 
 The idea is to generalize this code into Eigen. Hence it is written using Eigen.
 To use this function, download Eigen libraries from http://eigen.tuxfamily.org
 and install them locally. Use the headers:
 
 #include <Eigen/Sparse>
 #include <Eigen/IterativeLinearSolvers>
 #include <unsupported/Eigen/SparseExtra>
 #include <Eigen/Core>
 #include <Eigen/Dense>
 Use the bicgstabell() to solve the linear system of equations Ax=b. The arguments of the function is:
 1. A - mat
 2. b - rhs
 3. x - xguess
 4. Preconditioner (if any) - precond
 5. Number of iterations - iters
 6. Error - tol_error
 
 and compile the main file with g++ main.cpp -I /usr/local/include/eigen3.
