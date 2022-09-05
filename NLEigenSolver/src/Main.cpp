#include "nlpch.h"
#include "NLEigenSolver.h"

#include <stdio.h>

// Lapack wrappers
#include "f2c.h"
#include "clapack.h"

// Entry Point
int main(int argc, char* argv[])
{
    
    using Method = NLEigenMethods::Method;
    PROFILE_BEGIN_SESSION("Eigenvalue routine");
    std::shared_ptr<NLEigenSolver> app = NLEigenSolver::Create(Method::JacobiDavidson, argv[1]);
    bool status = app->execute();
    //bool status = app->findEigenvaluesFromInitialGuess();
    PROFILE_END_SESSION();

    // Check status
    if (!status)
    	return 1;


    // Checking the use of dstegr_
    // To do: 
     /*dstegr_(JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
         *ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
         *LIWORK, INFO);*/
    integer il = 0, iu = 0, nin = 6, nout = 6, info;
    char jobz = 'V';
    char range = 'A';
    double abstol = 1e-10;
    double vl = 0.0, vu = 0.0;
    integer m = 4, n = 4;
    integer ldz, liwork, lwork;
    ldz = n;
    liwork = 10 * n;
    lwork = 18 * n;
    double* diag = new double[n];
    double* subdiag = new double[n];
    double* eigenvalues = new double[n];
    double* work = new double[lwork];
    double* eigenvector = new double[ldz * n];
    integer* isuppz = new integer[(2 * n)];
    integer* iwork = new integer[liwork];
    
    // Set the values of the diagonal matrix
    /* 4x4 matrix T
     * 1 1 0 0
     * 1 4 2 0
     * 0 2 9 3
     * 0 0 3 16
     */
    diag[0] = 1.0;
    diag[1] = 4.0;
    diag[2] = -9.0;
    diag[3] = 16.0;
    subdiag[0] = 1.0;
    subdiag[1] = 2.0;
    subdiag[2] = 3.0;

    // Call DSTEGR subroutine (Algorithm to eigenlinear problem of symmetric tridiagonal  matrix)
    dstegr_(&jobz, &range, &n, diag, subdiag, &vl, &vu, &il, &iu, &abstol, &m, eigenvalues, eigenvector, &ldz, isuppz, work,
        &lwork, iwork, &liwork, &info);

    // Plot eigenvalues
    std::cout << "Eigenvalues: \n";
    std::cout << std::setprecision(16) << std::scientific;
   for (int i = 0; i < n; i++)
   {
        std::cout << eigenvalues[i] << " ";
   }
   std::cout << std::endl;

   std::cout << "Eigenvector: \n";
   for (int j = 0; j < n; j++)
   {
       for (int i = 0; i < ldz; i++)
       {
           std::cout << eigenvector[i*n+j] << " ";
       }
       std::cout << "\n";
   }

   // free memory
   delete[] diag;
   delete[] subdiag;
   delete[] eigenvalues;
   delete[] work;
   delete[] eigenvector;
   delete[] isuppz;
   delete[] iwork;


   return info;

}