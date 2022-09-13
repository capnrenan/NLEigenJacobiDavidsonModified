#pragma once
#include "NLEigenSolver.h"


class inverseFreeKrylovNLEigenSolver : public NLEigenSolver
{
	struct inverseKrylovData
	{
		int Dimensions;
		int NumberOfMassMtx;
		int NumberOfEigenvalues;
		int MaxIter;
		double Tolerance;
		//additional
		double Sigma;
	};

public:
	inverseFreeKrylovNLEigenSolver(const std::string& filepath);
	~inverseFreeKrylovNLEigenSolver();

	// Generalized inverse free Krylov method to compute the k smallest eigenpairs 
	// of the nonlinear symmetric eigenvalue problem
	virtual bool execute() override;

	//Miscellaneous
	virtual bool findEigenvaluesFromInitialGuess() override;

private:
	void readFileAndGetStiffMassMatrices(SparseMatrix& K0, std::vector<SparseMatrix>& MM, Vector& Omega);
	void printResults(Vector& Omega, DenseMatrix& Phi) const;
	void getFreqDependentStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, SparseMatrix& Kn, data_type omega);  // Kn(lr)
	void getFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mn, data_type omega);          
	void getGeneralizedFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mlrls, data_type lr, data_type ls);   // M(lr,ls)

	// Compute the smallest eigenpair using modified free inverse Krylov method
	void computeSmallestEigenpairKrylov(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM_Update,
										const std::vector<SparseMatrix>& MM_0, double& omega, Vector& phi, int numberOfKrylovBasis);

	// This functions compute the smallest eigenpairs of tridiagonal matrix 
	// by the MRRR algorithm
	void smallestEigenpairTriDiagMtx(double *diag, double* subdiag, int numberOfBasis, double& eigValue, Vector& eigVector);

	data_type RayleighQuotient(const SparseMatrix& A, const SparseMatrix& B, const Vector& eigVector);

	
	// Construct the Kyrlov subspace basis
	void orthoLanczosAlgorithm(double* diag, double* subdiag, const SparseMatrix& Ck, int numberOfBasis, const Vector& eigVector, DenseMatrix& Qm);
	void orthoArnoldiAlgorithm(const SparseMatrix& Ck, const SparseMatrix& B, int numberOfBasis, const Vector& eigVector, DenseMatrix& Zm);

	

private:
	std::string m_FilePath = "";
	inverseKrylovData* m_InputData = nullptr;
};


