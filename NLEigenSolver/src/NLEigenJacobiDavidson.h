#pragma once
#include "nlpch.h"

class NLEigenJacobiDavidson
{
	// Properties
public:
	unsigned int m_dimensions;
	unsigned int m_NumberOfMassMtx;
	unsigned int m_numberOfEigenValues;
	unsigned int m_MaxIter;
	double m_TOL;
private:

	//Methods
public:
	NLEigenJacobiDavidson(const char* filepath);
	~NLEigenJacobiDavidson();

	void execute();

private:
	void readData(const char* filepath);
	void printResults(Eigen::VectorXd& Omega, Eigen::MatrixXd& Phi) const;
    
	void getFreqDependentStiffMtx(const Eigen::MatrixXd& K0, const Eigen::MatrixXd& MM, Eigen::MatrixXd& Kn, double  omega);  // Kn(lr)
	void getFreqDependentMassMtx(const Eigen::MatrixXd& MM, Eigen::MatrixXd& Mn, double omega);                               // Mn(lr)
	void getGeneralizedFreqDependentMassMtx(const Eigen::MatrixXd& MM, Eigen::MatrixXd& Mlrls, double lr, double ls);         // M(lr,ls)
	void getEffectiveStiffMtx(const Eigen::MatrixXd& K0, const Eigen::MatrixXd& MM, Eigen::MatrixXd& Keff, double omega);     // Keff = K(lr)-lr*M(lr)
	void projectEffectiveStiffMatrix(Eigen::MatrixXd& Keff, Eigen::MatrixXd& B_s, unsigned int indexEig);
	bool iterativeLinearSolver(Eigen::MatrixXd& A, Eigen::VectorXd& b, Eigen::VectorXd& x);					// Solve Ax=b
};

