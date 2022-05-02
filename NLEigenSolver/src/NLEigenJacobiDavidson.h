#pragma once
#include "nlpch.h"

class NLEigenJacobiDavidson
{
	// Properties
private:
	int m_Dimensions;
	int m_NumberOfMassMtx;
	int m_NumberOfEigenValues;
	int m_MaxIter;
	double m_TOL;
	std::string m_FilePath;

public:
	NLEigenJacobiDavidson(const std::string& filepath);
	~NLEigenJacobiDavidson();

	bool execute();

private:
#ifdef ENABLE_BLAZE
	void readFileAndGetStiffMassMatrices(blaze::DynamicMatrix<double, blaze::rowMajor>& K0, std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM);
	void printResults(blaze::DynamicVector<double, blaze::columnVector>& Omega, blaze::DynamicMatrix<double, blaze::rowMajor>& Phi) const;
    
	void getFreqDependentStiffMtx(const blaze::DynamicMatrix<double, blaze::rowMajor>& K0, const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Kn, double omega);  // Kn(lr)
	void getFreqDependentMassMtx(const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Mn, double omega);                               // Mn(lr)
	void getGeneralizedFreqDependentMassMtx(const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Mlrls, double lr, double ls);         // M(lr,ls)
	void getEffectiveStiffMtx(const blaze::DynamicMatrix<double, blaze::rowMajor>& K0, const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Keff, double omega);     // Keff = K(lr)-lr*M(lr)
	void projectEffectiveStiffMatrix(blaze::DynamicMatrix<double, blaze::rowMajor>& Keff, blaze::DynamicMatrix<double, blaze::rowMajor>& B_s, int indexEig);
	bool iterativeLinearSolver(const blaze::DynamicMatrix<double, blaze::rowMajor>& A, const blaze::DynamicVector<double, blaze::columnVector>& b, blaze::DynamicVector<double, blaze::columnVector>& x);
#else
	void readFileAndGetStiffMassMatrices(Eigen::MatrixXd& K0, std::vector<Eigen::MatrixXd>& MM);
	void printResults(Eigen::VectorXd& Omega, Eigen::MatrixXd& Phi) const;

	void getFreqDependentStiffMtx(const Eigen::MatrixXd& K0, const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Kn, double omega);  // Kn(lr)
	void getFreqDependentMassMtx(const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Mn, double omega);                               // Mn(lr)
	void getGeneralizedFreqDependentMassMtx(const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Mlrls, double lr, double ls);         // M(lr,ls)
	void getEffectiveStiffMtx(const Eigen::MatrixXd& K0, const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Keff, double omega);     // Keff = K(lr)-lr*M(lr)
	void projectEffectiveStiffMatrix(Eigen::MatrixXd& Keff, Eigen::MatrixXd& B_s, int indexEig);
	bool iterativeLinearSolver(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);
#endif
	// Solve Ax=b
};

