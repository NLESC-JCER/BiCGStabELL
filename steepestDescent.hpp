//
//  steepestDescent.hpp
//  conjugateGradient
//
//  Created by Adithya Vijaykumar on 08/02/2019.
//  Copyright © 2019 Adithya Vijaykumar. All rights reserved.
//

#ifndef EIGEN_STEEPESTDESCENT
#define EIGEN_STEEPESTDESCENT

namespace Eigen {
    
    namespace internal {
        
        template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
        bool steepestdescent(const MatrixType& mat, const Rhs& rhs, Dest& x,
                             const Preconditioner& precond, Index& iters,
                             typename Dest::RealScalar& tol_error)
        {
            using std::sqrt;
            using std::abs;
            typedef typename Dest::RealScalar RealScalar;
            typedef typename Dest::Scalar Scalar;
            typedef Matrix<Scalar,Dynamic,1> VectorType;
            RealScalar tol = tol_error;
            Index maxIters = iters;
            
            Index n = mat.cols();
            
            VectorType r  = rhs - mat * x;
            VectorType r0 = r;
            
            RealScalar r0_sqnorm = r0.squaredNorm();
            RealScalar rhs_sqnorm = rhs.squaredNorm();
            if(rhs_sqnorm == 0)
            {
                x.setZero();
                return true;
            }
            
            
            RealScalar tol2 = tol*tol*rhs_sqnorm;
            Index i = 0;
            Index restarts = 0;
            
            while ( r.squaredNorm() > tol2 && i<iters )
            {

                VectorType r  = rhs - mat * x;
                std::cout << x(0) << " " << x(1) << std::endl;
                RealScalar alpha = r.dot(r)/(r.dot(mat*r));
                x = x + alpha*r;
                ++i;
                
            }
            tol_error = sqrt(r.squaredNorm()/rhs_sqnorm);
            iters = i;
            return true;
        }
        
    }
    
    template< typename _MatrixType,
    typename _Preconditioner = DiagonalPreconditioner<typename _MatrixType::Scalar> >
    class SteepestDescent;
    
    namespace internal {
        
        template< typename _MatrixType, typename _Preconditioner>
        struct traits<SteepestDescent<_MatrixType,_Preconditioner> >
        {
            typedef _MatrixType MatrixType;
            typedef _Preconditioner Preconditioner;
        };
        
    }
    
    template< typename _MatrixType, typename _Preconditioner>
    class SteepestDescent : public IterativeSolverBase<SteepestDescent<_MatrixType,_Preconditioner> >
    {
        typedef IterativeSolverBase<SteepestDescent> Base;
        using Base::matrix;
        using Base::m_error;
        using Base::m_iterations;
        using Base::m_info;
        using Base::m_isInitialized;
    public:
        typedef _MatrixType MatrixType;
        typedef typename MatrixType::Scalar Scalar;
        typedef typename MatrixType::RealScalar RealScalar;
        typedef _Preconditioner Preconditioner;
        
    public:
        
        /** Default constructor. */
        SteepestDescent() : Base() {}
        
        /** Initialize the solver with matrix \a A for further \c Ax=b solving.
         *
         * This constructor is a shortcut for the default constructor followed
         * by a call to compute().
         *
         * \warning this class stores a reference to the matrix A as well as some
         * precomputed values that depend on it. Therefore, if \a A is changed
         * this class becomes invalid. Call compute() to update it with the new
         * matrix A, or modify a copy of A.
         */
        template<typename MatrixDerived>
        explicit SteepestDescent(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}
        
        ~SteepestDescent() {}
        
        /** \internal */
        /** Loops over the number of columns of b and does the following:
         1. sets the tolerence and maxIterations
         2. Calls the function that has the core solver routine
         */
        template<typename Rhs,typename Dest>
        void _solve_with_guess_impl(const Rhs& b, Dest& x) const
        {
            bool failed = false;
            for(Index j=0; j<b.cols(); ++j)
            {
                //m_iterations = Base::maxIterations();
                //******************MANUALLY SET NUM ITERATIONS
                m_iterations = 30;
                m_error = Base::m_tolerance;
                
                typename Dest::ColXpr xj(x,j);
                if(!internal::steepestdescent(matrix(), b.col(j), xj, Base::m_preconditioner, m_iterations, m_error))
                    failed = true;
            }
            m_info = failed ? NumericalIssue
            : m_error <= Base::m_tolerance ? Success
            : NoConvergence;
            m_isInitialized = true;
        }
        
        /** \internal */
        /** Resizes the x vector to match the dimenstion of b and sets the elements to zero*/
        using Base::_solve_impl;
        template<typename Rhs,typename Dest>
        void _solve_impl(const MatrixBase<Rhs>& b, Dest& x) const
        {
            x.resize(this->rows(),b.cols());
            x.setZero();
            _solve_with_guess_impl(b,x);
        }
        
    protected:
        
    };
    
}
#endif /* steepestDescent_hpp */
