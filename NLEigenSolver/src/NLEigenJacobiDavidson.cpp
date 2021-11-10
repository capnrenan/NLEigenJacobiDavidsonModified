#include "nlpch.h"
#include "NLEigenJacobiDavidson.h"

#include <Eigen/IterativeLinearSolvers>

NLEigenJacobiDavidson::NLEigenJacobiDavidson(const char* filepath)
	: m_dimensions(0), m_NumberOfMassMtx(1), m_numberOfEigenValues(0),
	m_MaxIter(20), m_TOL(1e-12)
{

}

NLEigenJacobiDavidson::~NLEigenJacobiDavidson()
{

}

void NLEigenJacobiDavidson::execute()
{

}

void NLEigenJacobiDavidson::readData(const char* filepath)
{

}

void NLEigenJacobiDavidson::printResults(Eigen::VectorXd& Omega, Eigen::MatrixXd& Phi) const
{

}

void NLEigenJacobiDavidson::getFreqDependentStiffMtx(const Eigen::MatrixXd& K0, const Eigen::MatrixXd& MM, Eigen::MatrixXd& Kn, double omega)
{

}

void NLEigenJacobiDavidson::getFreqDependentMassMtx(const Eigen::MatrixXd& MM, Eigen::MatrixXd& Mn, double omega)
{

}

void NLEigenJacobiDavidson::getGeneralizedFreqDependentMassMtx(const Eigen::MatrixXd& MM, Eigen::MatrixXd& Mlrls, double lr, double ls)
{

}

void NLEigenJacobiDavidson::getEffectiveStiffMtx(const Eigen::MatrixXd& K0, const Eigen::MatrixXd& MM, Eigen::MatrixXd& Keff, double omega)
{

}

void NLEigenJacobiDavidson::projectEffectiveStiffMatrix(Eigen::MatrixXd& Keff, Eigen::MatrixXd& B_s, unsigned int indexEig)
{

}

bool NLEigenJacobiDavidson::iterativeLinearSolver(Eigen::MatrixXd& A, Eigen::VectorXd& b, Eigen::VectorXd& x)
{
	// Set the iterative linear solver (Conjugate Gradients)
	// The number of max. of iter
	Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower | Eigen::Upper> linsolver;
	linsolver.setTolerance(1e-12);
	linsolver.compute(A);

	//Solve
	x = linsolver.solve(b);
	
	//LOG_INFO("#Iteration: {0}   Estimated error: {1}", linsolver.iterations(), linsolver.error());

	bool status = true;
	
	if (linsolver.error() < 1e-12)
		status = false;
	
	return status;
}

