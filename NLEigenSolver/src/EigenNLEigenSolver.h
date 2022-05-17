#pragma once
#include "NLEigenSolver.h"

#include <Eigen/Dense>
#include <Eigen/Core>

class EigenNLEigenSolver : public NLEigenSolver
{
public:
	EigenNLEigenSolver(const std::string& filepath);
	~EigenNLEigenSolver();

	virtual bool execute() override;

private:
	void readFileAndGetStiffMassMatrices(Eigen::MatrixXd& K0, std::vector<Eigen::MatrixXd>& MM);
	void printResults(Eigen::VectorXd& Omega, Eigen::MatrixXd& Phi) const;
	void getFreqDependentStiffMtx(const Eigen::MatrixXd& K0, const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Kn, double omega);  // Kn(lr)
	void getFreqDependentMassMtx(const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Mn, double omega);                               // Mn(lr)
	void getGeneralizedFreqDependentMassMtx(const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Mlrls, double lr, double ls);         // M(lr,ls)
	void getEffectiveStiffMtx(const Eigen::MatrixXd& K0, const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Keff, double omega);     // Keff = K(lr)-lr*M(lr)
	void projectEffectiveStiffMatrix(Eigen::MatrixXd& Keff, Eigen::MatrixXd& B_s, int indexEig);
	bool iterativeLinearSolver(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);
	void UpdateEigenvectorSolution(Eigen::MatrixXd& Keff, Eigen::MatrixXd& Phi, Eigen::MatrixXd& B_r, int index);

public:
	int m_Dimensions;
	int m_NumberOfMassMtx;
	int m_NumberOfEigenValues;
	int m_MaxIter;
	double m_TOL;
	std::string m_FilePath;
};