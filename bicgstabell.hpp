//
//  bicgstabell.h
//  conjugateGradient
//
//  Created by Adithya Vijaykumar on 19/02/2019.
//  Copyright © 2019 Adithya Vijaykumar. All rights reserved.
//

#ifndef bicgstabell_h
#define bicgstabell_h

#include <vector>

namespace Eigen {
    
    namespace internal {
        
        template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
        bool bicgstabell(const MatrixType& mat, const Rhs& rhs, Dest& xguess,
                                const Preconditioner& precond, Index& iters,
                         typename Dest::RealScalar& tol_error)
        {
            
            using std::sqrt;
            using std::abs;
            typedef typename Dest::RealScalar RealScalar;
            typedef typename Dest::Scalar Scalar;
            int l=5;
            //start with k=-l or k+l=0
            //TODO Set l=2 as a default**************************
            int L = l;
            int k = -L;
            
            RealScalar tol = tol_error;
            Index maxIters = iters;

            typedef Matrix<Scalar, Dynamic,1> VectorType;
            typedef Matrix<RealScalar, Dynamic, Dynamic> DenseMatrixType;
           
            //We start with an initial guess x_0 and let us set r_0 as (residual calculated from x_0)
            VectorType x0 = xguess;
            VectorType r0  = rhs - mat * x0; //r_0
            
            VectorType u0;
            u0.setZero();
            VectorType rShadow = r0; //shadow of r0 is r0
            
            //write over these variable the successive iteratively got values
            VectorType x = x0;
            VectorType u = u0;
            
            RealScalar deltaNew = r0.squaredNorm();//r.r
            RealScalar delta0 = deltaNew;
            RealScalar rhs_sqnorm = rhs.squaredNorm();
            
            //Other vectors and scalars initialisation
            
            RealScalar rho0 = 1.0;
            RealScalar alpha = 0.0;
            RealScalar omega = 1.0;
            
            
            //Consequence of A being invertible if Ax=0 ==> x=0
            if(rhs_sqnorm == 0)
            {
                x.setZero();
                return true;
            }
            
            
            RealScalar tol2 = tol*tol*delta0;
            Index currIter = 0;
            
            std::vector<VectorType> rHat(L+1); //rj hat is a vector of vectors!!!
            std::vector<VectorType> uHat(L+1);
            //std::vector<VectorType> xHat(L);
            rHat[0] = r0;
            uHat[0].setZero(r0.rows());
            
            while ( rHat[0].squaredNorm() > tol2 && currIter<maxIters )
            {
                k = k+L;
                rho0 = -omega*rho0;
                //std::cout << " uhat0:" <<uHat[0] << " rhat[0] " << rHat[0] << "  x " << x;
                //std::cout << " rho0: " << rho0 << " alpha: " << alpha << std::endl;
                for(Index j=0;j<=L-1;++j)   //THIS FOR LOOP IS THE BI-CG PART
                {
                    //std::cout <<"j:" << j << "\n rhat: \n" << rHat[j] << std::endl;
                    RealScalar rho1 = rHat[j].dot(rShadow);
                    //std::cout << "rho1: " << rho1 << std::endl;
                    RealScalar beta = alpha * (rho1/rho0);
                    rho0 = rho1;
                    for(Index i=0; i<=j; ++i)
                    {
                        
                        uHat[i] = rHat[i] - beta*uHat[i];
                        
                    }
                    //std::cout<<"uhatj"<<uHat[j]<<std::endl;
                    uHat[j+1] = mat * uHat[j];
                    //std::cout << "alpha: " << alpha << std::endl;
                    alpha = rho0/(uHat[j+1].dot(rShadow));
                    //std::cout<<"uhatj_afterA"<<uHat[j+1]<<std::endl;
                    
                    for(Index i=0; i<=j; ++i)
                    {
                        rHat[i] = rHat[i] - alpha*uHat[i+1];
                    }
                    //std::cout << "rHat[0]: " << rHat[0] << std::endl;
                    rHat[j+1] = mat * rHat[j];
                   //std::cout << "rHat[j+1]: " << rHat[j+1] << std::endl;
                    x = x + alpha*uHat[0];
                    
                }
                //std::cout<<"x after bcig:"<<x<<std::endl;
                //std::cout<<"alpha:"<<alpha<<std::endl;
                //std::cout<<"uhat0:"<<uHat[0]<<std::endl;
                
                //THE MINIMAL RESIDUAL PART STARTS NOW
                
                DenseMatrixType tau(L, L+1);
                std::vector<RealScalar> sigma(L+1);
                std::vector<RealScalar> gammaP(L+1);
                
                for(Index j=1;j<=L;++j)
                {
                    for(Index i=1;i<=j-1;++i)
                    {
                        
                        tau(i,j) = rHat[j].dot(rHat[i])/sigma[i];
                        rHat[j] = rHat[j] - tau(i,j) * rHat[i];
                        
                    }
                    
                    sigma[j] = rHat[j].dot(rHat[j]);
                    gammaP[j] = rHat[0].dot(rHat[j])/sigma[j];
                    
                }
                
                std::vector<RealScalar> gamma(L+1);
                gamma[L] = gammaP[L];
                omega = gamma[L];
                
                for(Index j=L-1;j>=1; --j)
                {
                    RealScalar sum = 0.0;
                    for(Index i=j+1; i<=L; ++i)
                    {
                        sum += tau(j,i)*gamma[i];
                    }
                    gamma[j] = gammaP[j] - sum;
                }
                
                std::vector<RealScalar> gammaPP(L);
                for(Index j=1; j<=L-1; ++j)
                {
                    RealScalar sum = 0.0;
                    for(Index i=j+1; i<=L-1; ++i)
                    {
                        sum += tau(j,i)*gamma[i+1];
                    }
                    gammaPP[j] = gamma[j+1] + sum;
                }
                
                x = x + gamma[1] * rHat[0];
                rHat[0] = rHat[0] - gammaP[L]*rHat[L];
                uHat[0] = uHat[0] - gamma[L]*uHat[L];
                
                //std::cout<<"x at end:"<<x<<std::endl;
                
                for (Index j=1; j<=L-1; ++j)
                {
                    uHat[0] = uHat[0] - gamma[j]*uHat[j];
                    x = x + gammaPP[j]*rHat[j];
                    rHat[0] = rHat[0] - gammaP[j]*rHat[j];
                    
                    
                }
                ++currIter;
                /*std::cout<<"gammaPP"<<gammaPP[1]<<std::endl;
                std::cout<<"x after loop:"<<x<<std::endl;
                std::cout<<"tau:"<<tau<<std::endl;
                std::cout<<"gammaP"<<gammaP[0]<<" "<<gammaP[1]<<" "<<gammaP[2]<<std::endl;
                std::cout<<"gammaPP"<<gammaPP[0]<<" "<<gammaPP[1]<<" "<<gammaPP[2]<<std::endl;
                std::cout<<"gammas"<<gamma[0]<<" "<<gamma[1]<<" "<<gamma[2]<<std::endl;
                *///return 0;
            }
            
            tol_error = sqrt(rHat[0].squaredNorm()/rhs_sqnorm);
            iters = currIter;
            xguess = x;
            return true;
        }
        
    }
    
    template< typename _MatrixType,
    typename _Preconditioner = DiagonalPreconditioner<typename _MatrixType::Scalar> >
    class BicgstabEll;
    
    namespace internal {
        
        template< typename _MatrixType, typename _Preconditioner>
        struct traits<BicgstabEll<_MatrixType,_Preconditioner> >
        {
            typedef _MatrixType MatrixType;
            typedef _Preconditioner Preconditioner;
        };
        
    }
    
    template< typename _MatrixType, typename _Preconditioner>
    class BicgstabEll : public IterativeSolverBase<BicgstabEll<_MatrixType,_Preconditioner> >
    {
        typedef IterativeSolverBase<BicgstabEll> Base;
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
        BicgstabEll() : Base() {}
        
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
        explicit BicgstabEll(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}
        
        ~BicgstabEll() {}
        
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
                m_iterations = Base::maxIterations();
                //******************MANUALLY SET NUM ITERATIONS
                //m_iterations = 30;
                m_error = Base::m_tolerance;
                
                typename Dest::ColXpr xj(x,j);
                if(!internal::bicgstabell(matrix(), b.col(j), xj, Base::m_preconditioner, m_iterations, m_error))
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

#endif /* bicgstabell_h */

