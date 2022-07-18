#pragma once
#include "NLEigenSolver.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <unsupported/Eigen/MPRealSupport>



#define QUAD_PRECISION 0
// Check the quad precision
#if QUAD_PRECISION
	using data_type = mpfr::mpreal;
	using DenseMatrix = Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic>;
	using EigenVector = Eigen::Vector<mpfr::mpreal, Eigen::Dynamic>;

	using SparseMatrix = Eigen::SparseMatrix<mpfr::mpreal>;

#else
	using data_type = double;
	using DenseMatrix = Eigen::MatrixXd;
	using EigenVector = Eigen::VectorXd;

	using SparseMatrix = Eigen::SparseMatrix<data_type>;
	
#endif

// Old implementation
#if OLD_IMPL
class EigenNLEigenSolver : public NLEigenSolver
{
public:
	EigenNLEigenSolver(const std::string& filepath);
	~EigenNLEigenSolver();

	virtual bool execute() override;

private:
	void readFileAndGetStiffMassMatrices(DenseMatrix& K0, std::vector<DenseMatrix>& MM, EigenVector& Omega);
	void printResults(EigenVector& Omega, DenseMatrix& Phi) const;
	void getFreqDependentStiffMtx(const DenseMatrix& K0, const std::vector<DenseMatrix>& MM, DenseMatrix& Kn, data_type omega);  // Kn(lr)
	void getFreqDependentMassMtx(const std::vector<DenseMatrix>& MM, DenseMatrix& Mn, data_type omega);                               // Mn(lr)
	void getGeneralizedFreqDependentMassMtx(const std::vector<DenseMatrix>& MM, DenseMatrix& Mlrls, data_type lr, data_type ls);         // M(lr,ls)
	void getEffectiveStiffMtx(const DenseMatrix& K0, const std::vector<DenseMatrix>& MM, DenseMatrix& Keff, data_type omega);     // Keff = K(lr)-lr*M(lr)
	void projectEffectiveStiffMatrix(DenseMatrix& Keff, DenseMatrix& B_s, int indexEig);
	bool iterativeLinearSolver(const DenseMatrix& A, const EigenVector& b, EigenVector& x);
	void UpdateEigenvectorSolution(DenseMatrix& Keff, DenseMatrix& Phi, DenseMatrix& B_r, int index);

private:
	int m_Dimensions;
	int m_NumberOfMassMtx;
	int m_NumberOfEigenValues;
	int m_MaxIter;
	double m_TOL;
	std::string m_FilePath;
	bool m_hasInitialTrial = false;
};
#else
	//  Implementation with support to sparse matrix (Optmizing the code!!!!)
	class EigenNLEigenSolver : public NLEigenSolver
	{
	public:
		EigenNLEigenSolver(const std::string& filepath);
		~EigenNLEigenSolver();

		virtual bool execute() override;
		//Miscellaneous
		virtual bool findEigenvaluesFromInitialGuess() override;

	private:
		void readFileAndGetStiffMassMatrices(SparseMatrix& K0, std::vector<SparseMatrix>& MM, EigenVector& Omega);
		void printResults(EigenVector& Omega, DenseMatrix& Phi) const;
		void getFreqDependentStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, SparseMatrix& Kn, data_type omega);  // Kn(lr)
		void getFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mn, data_type omega);                               // Mn(lr)
		void getGeneralizedFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mlrls, data_type lr, data_type ls);         // M(lr,ls)
		void getEffectiveStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, DenseMatrix& Keff, data_type omega);     // Keff = K(lr)-lr*M(lr)
		void projectEffectiveStiffMatrix(DenseMatrix& Keff, DenseMatrix& B_s, int indexEig);
		bool iterativeLinearSolver(const DenseMatrix& A, const EigenVector& b, EigenVector& x);
		void UpdateEigenvectorSolution(DenseMatrix& Keff, DenseMatrix& Phi, DenseMatrix& B_r, int index);


	private:
		int m_Dimensions;
		int m_NumberOfMassMtx;
		int m_NumberOfEigenValues;
		int m_MaxIter;
		double m_TOL;
		std::string m_FilePath;
		bool m_hasInitialTrial = false;

	};
#endif

