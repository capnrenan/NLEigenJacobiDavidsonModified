#include "nlpch.h"
#include "NLEigenJacobiDavidson.h"

#include <Eigen/IterativeLinearSolvers>

NLEigenJacobiDavidson::NLEigenJacobiDavidson(const std::string& filepath)
	: m_Dimensions(0), m_NumberOfMassMtx(1), m_NumberOfEigenValues(0),
	m_MaxIter(20), m_TOL(1e-12), m_FilePath(filepath)
{
	// Initialize the log system
	Log::Init();
}

NLEigenJacobiDavidson::~NLEigenJacobiDavidson()
{

}

void NLEigenJacobiDavidson::execute()
{
	//Reading data
	LOG_INFO("Reading filedata...\n");
	Eigen::MatrixXd K0;
	std::vector<Eigen::MatrixXd> MM;
	readFileAndGetStiffMassMatrices(K0, MM);


}

void NLEigenJacobiDavidson::readFileAndGetStiffMassMatrices(Eigen::MatrixXd& K0, std::vector<Eigen::MatrixXd>& MM)
{

}

void NLEigenJacobiDavidson::printResults(Eigen::VectorXd& Omega, Eigen::MatrixXd& Phi) const
{

}

void NLEigenJacobiDavidson::getFreqDependentStiffMtx(const Eigen::MatrixXd& K0, const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Kn, double omega)
{
	//Initialize
	Kn = K0;
    
	for (int jj = 1; jj < m_NumberOfMassMtx; jj++)
	{

		Kn += (jj)*pow(omega, jj + 1.0) * MM[jj];
	}
}

void NLEigenJacobiDavidson::getFreqDependentMassMtx(const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Mn, double omega)
{
	//Initialize
	Mn = MM[0];

	for (int jj = 1; jj < m_NumberOfMassMtx; jj++)
	{

		Mn += (jj + 1.0) * pow(omega, jj) * MM[jj];
	}
	
}

void NLEigenJacobiDavidson::getGeneralizedFreqDependentMassMtx(const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Mlrls, double lr, double ls)
{
	//Initialize
	Mlrls.setZero();

	for (int jj = 0; jj < m_NumberOfMassMtx; jj++)
	{
		for (int kk = 0; kk < jj + 1; kk++)
		{
			Mlrls += pow(lr, kk) * pow(ls, jj - kk) * MM[jj];
		}
	}
}

void NLEigenJacobiDavidson::getEffectiveStiffMtx(const Eigen::MatrixXd& K0, const std::vector<Eigen::MatrixXd>& MM, Eigen::MatrixXd& Keff, double omega)
{
	// Initialize
	Keff = K0;

	for (int jj = 0; jj < m_NumberOfMassMtx; jj++)
	{
		Keff -= pow(omega, jj + 1.0) * MM[jj];
	}
}

void NLEigenJacobiDavidson::projectEffectiveStiffMatrix(Eigen::MatrixXd& Keff, Eigen::MatrixXd& B_s, int indexEig)
{
	// Project the effective stiffness matrix onto the subspace
	// orthogonal to all preceding eigenvectors and add the orthogonal
	// projector to make it nonsingular
	for (int ii = 0; ii < indexEig; ii++)
	{
		Keff += (B_s.col(ii) - Keff * B_s.col(ii)) * (B_s.col(ii).transpose());
	}
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
	{
		status = false;
		LOG_ASSERT(status,"The iterative linear solver has reach the max. iterations with error below of the tolerance!")
	}
		
	return status;
}

