#pragma once
#include "NLEigenSolver.h"

#include <Eigen/Dense>
#include <Eigen/Core>

#define QUAD_PRECISION 1

// Check the quad precision
#if QUAD_PRECISION
	using EigenMatrix = Eigen::MatrixX<long double>;
	using EigenVector = Eigen::VectorX<long double>;

	using data_type = long double;

#else
	using EigenMatrix = Eigen::MatrixXd;
	using EigenVector = Eigen::VectorXd;

	using data_type =  double;
#endif


class EigenNLEigenSolver : public NLEigenSolver
{
public:
	EigenNLEigenSolver(const std::string& filepath);
	~EigenNLEigenSolver();

	virtual bool execute() override;

private:
	void readFileAndGetStiffMassMatrices(EigenMatrix& K0, std::vector<EigenMatrix>& MM, EigenVector& Omega);
	void printResults(EigenVector& Omega, EigenMatrix& Phi) const;
	void getFreqDependentStiffMtx(const EigenMatrix& K0, const std::vector<EigenMatrix>& MM, EigenMatrix& Kn, data_type omega);  // Kn(lr)
	void getFreqDependentMassMtx(const std::vector<EigenMatrix>& MM, EigenMatrix& Mn, data_type omega);                               // Mn(lr)
	void getGeneralizedFreqDependentMassMtx(const std::vector<EigenMatrix>& MM, EigenMatrix& Mlrls, data_type lr, data_type ls);         // M(lr,ls)
	void getEffectiveStiffMtx(const EigenMatrix& K0, const std::vector<EigenMatrix>& MM, EigenMatrix& Keff, data_type omega);     // Keff = K(lr)-lr*M(lr)
	void projectEffectiveStiffMatrix(EigenMatrix& Keff, EigenMatrix& B_s, int indexEig);
	bool iterativeLinearSolver(const EigenMatrix& A, const EigenVector& b, EigenVector& x);
	void UpdateEigenvectorSolution(EigenMatrix& Keff, EigenMatrix& Phi, EigenMatrix& B_r, int index);

private:
	int m_Dimensions;
	int m_NumberOfMassMtx;
	int m_NumberOfEigenValues;
	int m_MaxIter;
	double m_TOL;
	std::string m_FilePath;
	bool m_hasInitialTrial = false;
};