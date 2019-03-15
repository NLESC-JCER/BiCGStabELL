//
//  main.cpp
//  conjugateGradient
//
//  Created by Adithya Vijaykumar on 08/02/2019.
//  Copyright © 2019 Adithya Vijaykumar. All rights reserved.
//


#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "steepestDescent.hpp"
#include "conjugateGradient.hpp"
#include "bicgstabell.hpp"
using namespace std ;
using namespace Eigen ;
typedef Eigen::SparseMatrix<double> SpMat;
int main ()
{
    SpMat A;
    VectorXd b,x,t,guess(2,1);
    guess(0) = -2;
    guess(1) = -2;
    //loadMarket(A, "data/3x3/A_3x3.txt");
    //loadMarketVector(b,"data/3x3/b_3x3.txt");
    
    loadMarket(A, "test_matrices/sleijpen_example_A.txt");
    loadMarketVector(b,"test_matrices/sleijpen_example_b.txt");
    //std::cout << A << std::endl;
    //std::cout << b << std::endl;
    
    BicgstabEll<SparseMatrix<double> > solver;
    solver.compute(A);
    x=solver.solve(b);
    //x=solver.solveWithGuess(b, guess);
    //std::cout << "#iterations:     " << solver.iterations() << std::endl;
    //std::cout << "estimated error: " << solver.error()      << std::endl;
    /* ... update b ... */
    //x = solver.solve(b); // solve again
    /*std::cout << A << std::endl;
     std::cout << b << std::endl;
     std::cout << x<< std::endl;*/
    /*std::cout << x<< std::endl;
     BiCGSTAB<SparseMatrix<double> > solver;
     solver.compute(A);
     x=solver.solve(b);*/
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
    /* ... update b ... */
    //x = solver.solve(b); // solve again
    //std::cout << A << std::endl;
    //std::cout << b << std::endl;
    std::cout << x << std::endl;
    return 0;
}
