#pragma once
#include "NLEigenSolver.h"


class inverseFreeKrylovNLEigenSolver : public NLEigenSolver
{
	struct resultStatus
	{
		data_type Convergence;
		int NumberOfIterations;
		bool Status; // true if the result converge
	};

	// It is only used the modified Incomplete Cholesky with
	// dual threshold to compute the precondioner as available
	// in Eigen 3.4.90. 
	// 
	// Ref.
	// 1. C-J. Lin and J. J. Moré, Incomplete Cholesky Factorizations with Limited memory, 
	// SIAM J. Sci. Comput. 21(1), pp. 24-45, 1999
	// 2. https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteCholesky.html#ad822e66656638a9cf84087ed228514e7
	// 
	struct preconditionerOptions
	{
		bool IsUsingPreconditioner = false;
		double InitialShiftSigma = 1e-3;
	};

	struct inverseKrylovData
	{
		int Dimensions;
		int NumberOfMassMtx;
		int NumberOfEigenvalues = 1;
		int MaxIter = 100;
		double Tolerance = 1e-10;
		//additional
		double Sigma = 1e6;
		int NumberOfKrylovBasis = 10;

		preconditionerOptions PrecondOptions;
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
	void readFileAndGetStiffMassMatrices(SparseMatrix& K0, std::vector<SparseMatrix>& MM);
	void printInfo();
	void printResults(Vector& Omega, DenseMatrix& Phi) const;
	
	template<typename MtxType>
	void getFreqDependentStiffMtx(const MtxType& K0, const std::vector<MtxType>& MM, MtxType& Kn, data_type omega);  // Kn(lr)
	
	template<typename MtxType>
	void getFreqDependentMassMtx(const std::vector<MtxType>& MM, MtxType& Mn, data_type omega);

	template<typename MtxType>
	void getGeneralizedFreqDependentMassMtx(const std::vector<MtxType>& MM, MtxType& Mlrls, data_type lr, data_type ls);   // M(lr,ls)

	// Generalized deflation technique for nonlinear eigenvalue problems
	void generalizedDeflationProcedure(DenseMatrix& K0, std::vector<DenseMatrix>& MM_Update, DenseMatrix& Mlrls, data_type omega, Vector& phi);

	// Compute the smallest eigenpair using modified free inverse Krylov method
	void computeSmallestEigenpairKrylov(const DenseMatrix& K0, const std::vector<DenseMatrix>& MM_Update,
										const std::vector<SparseMatrix>& MM_0, data_type& omega, Vector& phi,
										int numberOfKrylovBasis, resultStatus& status);

	// This functions compute the smallest eigenpairs of tridiagonal matrix 
	// by the MRRR algorithm
	void smallestEigenpairTriDiagMtx(double *diag, double* subdiag, int numberOfBasis, double& eigValue, Vector& eigVector);

	template<typename MtxType>
	data_type RayleighQuotient(const MtxType& A, const MtxType& B, const Vector& eigVector);

	// Construct the Kyrlov subspace basis
	void orthoLanczosAlgorithm(double* diag, double* subdiag, const DenseMatrix& Ck, int numberOfBasis, const Vector& eigVector, DenseMatrix& Qm);
	void orthoArnoldiAlgorithm(const DenseMatrix& Ck, const DenseMatrix& B, int numberOfBasis, const Vector& eigVector, DenseMatrix& Zm);

	// preconditioned version of Arnoldi algorithm
	// invLktLk - inv(L^t*L)
	// L is lower triangle matrix obtained by a incomplemte Cholesky factorization. 
	void orthoPreconditionedArnoldiAlgorithm(const DenseMatrix& Ck, const DenseMatrix& B, int numberOfBasis, const Vector& eigVector, DenseMatrix& Zm, DenseMatrix& invLktLk);

	bool SetAndComputePreconditioner(const DenseMatrix& Ck, DenseMatrix& Minv);

private:
	std::string m_FilePath = "";
	inverseKrylovData m_InputData;
};


