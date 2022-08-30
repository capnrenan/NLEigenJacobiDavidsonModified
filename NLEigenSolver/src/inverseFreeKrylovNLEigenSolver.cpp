#include "nlpch.h"
#include "inverseFreeKrylovNLEigenSolver.h"

// Lapack wrappers
#include "f2c.h"
#include "clapack.h"

inverseFreeKrylovNLEigenSolver::inverseFreeKrylovNLEigenSolver(const std::string& filepath)
    : m_FilePath(filepath)
{


}

inverseFreeKrylovNLEigenSolver::~inverseFreeKrylovNLEigenSolver()
{
}

bool inverseFreeKrylovNLEigenSolver::execute()
{
	return false;
}

bool inverseFreeKrylovNLEigenSolver::findEigenvaluesFromInitialGuess()
{
	return false;
}

void inverseFreeKrylovNLEigenSolver::readFileAndGetStiffMassMatrices(SparseMatrix& K0, std::vector<SparseMatrix>& MM, Vector& Omega)
{
}

void inverseFreeKrylovNLEigenSolver::printResults(Vector& Omega, DenseMatrix& Phi) const
{
}

void inverseFreeKrylovNLEigenSolver::getFreqDependentStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, SparseMatrix& Kn, data_type omega)
{
}

void inverseFreeKrylovNLEigenSolver::getFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mn, data_type omega)
{
}

void inverseFreeKrylovNLEigenSolver::getGeneralizedFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mlrls, data_type lr, data_type ls)
{
}

void inverseFreeKrylovNLEigenSolver::smallestEigenpairTriDiagMtx(double* diag, double* subdiag, int numberOfBasis, double& eigValue, Vector& eigVector)
{
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
    integer m = numberOfBasis, n = numberOfBasis;
    integer ldz, liwork, lwork;
    ldz = n;
    liwork = 10 * n;
    lwork = 18 * n;
    double* eigenvalues = new double[n];
    double* work = new double[lwork];
    double* eigenvector = new double[ldz * n];
    integer* isuppz = new integer[(2 * n)];
    integer* iwork = new integer[liwork];

    // Call DSTEGR subroutine (Algorithm to eigenlinear problem of symmetric tridiagonal  matrix)
    // MRRR algorithm
    dstegr_(&jobz, &range, &n, diag, subdiag, &vl, &vu, &il, &iu, &abstol, &m, eigenvalues, eigenvector, &ldz, isuppz, work,
        &lwork, iwork, &liwork, &info);

    // Get the smallest eigenpair
    eigValue = eigenvalues[0];
    for (uint32_t i = 0; i < n; i++)
    {
        eigVector(i) = eigenvector[i];
    }

    // free memory
    delete[] eigenvalues;
    delete[] work;
    delete[] eigenvector;
    delete[] isuppz;
    delete[] iwork;

}

data_type inverseFreeKrylovNLEigenSolver::RayleighQuotient(const SparseMatrix& A, const SparseMatrix& B, const Vector& eigVector)
{
    data_type result = 0.0;
    data_type num = eigVector.transpose() * A * eigVector;
    data_type den = eigVector.transpose() * B * eigVector;
    result = num / den;
	return result;
}

void inverseFreeKrylovNLEigenSolver::orthoLanczosAlgorithm(double* diag, double* subdiag, const SparseMatrix& Ck, int numberOfBasis, const Vector& eigVector, DenseMatrix& Qm)
{

}

void inverseFreeKrylovNLEigenSolver::orthoArnoldiAlgorithm(const SparseMatrix& Ck, const SparseMatrix& B, int numberOfBasis, const Vector& eigVector, DenseMatrix& Zm)
{
}
