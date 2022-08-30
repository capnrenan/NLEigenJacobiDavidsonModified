#include "nlpch.h"
#include "MathUtils.h"
#include <Eigen/LU>
#include "JacobiDavidsonNLEigenSolver.h"

JacobiDavidsonNLEigenSolver::JacobiDavidsonNLEigenSolver(const std::string& filepath)
	: m_Dimensions(0), m_NumberOfMassMtx(1), m_NumberOfEigenValues(0),
	m_MaxIter(200), m_TOL(1e-14), m_FilePath(filepath)
{
	Eigen::initParallel();

	#if QUAD_PRECISION
		LOG_WARN("Quad precision is enabled!");
	#else
		LOG_WARN("Quad precision is unabled!");
	#endif


}



JacobiDavidsonNLEigenSolver::~JacobiDavidsonNLEigenSolver()
{

}

bool JacobiDavidsonNLEigenSolver::findEigenvaluesFromInitialGuess()
{
	LOG_INFO("Finding eigenvalues from an initial guess using Regula-Falsi\n");
	//Initializing and reading data
	LOG_INFO("Reading filedata...\n");
	SparseMatrix K0;
	std::vector<SparseMatrix> MM;
	Vector Omega;
	readFileAndGetStiffMassMatrices(K0, MM, Omega);
	bool status = true;
	DenseMatrix Ktemp(m_Dimensions, m_Dimensions);

	// Loop to find eigenvalue from a initial guess
	// using the linear regula-falsi algorithm
	for (int ie = 0; ie < m_NumberOfEigenValues; ie++)
	{
		LOG_INFO("---------------------------------\nEigenvalue #{0}:", ie);
		Ktemp.setZero();
		data_type value = MathUtils::linearRegularFalsi<data_type>(Omega(ie), [&](data_type value)
			{
				//getEffectiveStiffMtx(K0, MM, Ktemp, value);
				Ktemp = K0;
				for (int jj = 0; jj < m_NumberOfMassMtx; jj++)
				{
					Ktemp -= pow(value, jj + 1.0) * MM[jj];
				}
				auto sqrtDet = Ktemp.llt().matrixL().determinant();
				return sqrtDet * sqrtDet;
			}, status, m_TOL, m_MaxIter);

		// Update the eigenvalue
		Omega(ie) = value;
	}

	// Save the eigenproblem results
	LOG_INFO("Save the eigenvalues!");
	// Get the directory path
	std::string directory, resultFile;
	const size_t last_slash_idx = m_FilePath.rfind('\\');
	if (std::string::npos != last_slash_idx)
		directory = m_FilePath.substr(0, last_slash_idx);

	resultFile = directory + "\\Omega_Root.dat";
	std::ofstream out1;
	out1.open(resultFile);

	if (!out1)
		LOG_ASSERT(false, "ERROR: Error in opening the file!");

	//Save Omega
	out1 << std::setprecision(16) << std::scientific << Omega;
	out1.close();

	return status;
}

bool JacobiDavidsonNLEigenSolver::execute()
{
	//Reading data
	LOG_INFO("Reading filedata...\n");

	SparseMatrix K0;
	std::vector<SparseMatrix> MM;
	Vector Omega;
	readFileAndGetStiffMassMatrices(K0, MM, Omega);

	LOG_INFO("Initializing the matrices...\n");
	//Initialize the matrices and set them as zero
	DenseMatrix B_r(m_Dimensions, m_NumberOfEigenValues);
	DenseMatrix Keff(m_Dimensions, m_Dimensions);
	DenseMatrix Phi(m_Dimensions, m_NumberOfEigenValues);
	SparseMatrix Kn(m_Dimensions, m_Dimensions);
	SparseMatrix Mn(m_Dimensions, m_Dimensions);
	SparseMatrix Mlrls(m_Dimensions, m_Dimensions);

	LOG_INFO("Solving with TOL = {0}\n", m_TOL);

	//Set as zero
	B_r.setZero(); Keff.setZero(); Phi.setOnes();

	//Auxiliary variables
	data_type conv, normBr, PtMP, PtKP, theta;
	int iterK;

	LOG_INFO("Processing...");
	// Loop in each eigenvalue
	for (int ie = 0; ie < m_NumberOfEigenValues; ie++)
	{
		// Set the convergence parameters
		conv = 1.0;
		iterK = 0;

		LOG_INFO("---------------------------------\nEigenvalue #{0}:", ie);

		if (!m_hasInitialTrial)
		{
			if (ie > 0)
			{
				Omega(ie) = Omega(ie - 1);
			}
		}

		while (abs(conv) > m_TOL)
		{
			// Orthogonalized phi_r with respect to phi_s
			for (int is = 0; is < ie + 1; is++)
			{
				// Get the generalized frequency-dependent mass matrix
				getGeneralizedFreqDependentMassMtx(MM, Mlrls, Omega(ie), Omega(is));
				B_r.col(ie) = Mlrls * Phi.col(is);

				if (is > 0)
				{
					for (int el = 0; el < is; el++)
					{
						//Project b_s = b_s - b_el*(b_el.t()*b_s)
						data_type temp = (B_r.col(el).transpose() * B_r.col(is));
						B_r.col(is) += -B_r.col(el) * temp;
					}
				}

				// Normalize
				normBr = sqrt(B_r.col(is).transpose() * B_r.col(is));
				B_r.col(is) = 1.0 / (normBr)*B_r.col(is);
			}


			// Orthogonalize phi_e with respect to the preceding Vector phi
			if (ie > 0)
			{
				for (int is = 0; is < ie; is++)
				{
					Phi.col(ie) += -B_r.col(is) * (B_r.col(is).transpose() * Phi.col(ie));
				}
			}

			// Evaluate the effective stiffness matrix
			getEffectiveStiffMtx(K0, MM, Keff, Omega(ie));
			

			// Update Vector solution
			std::thread thread0(&JacobiDavidsonNLEigenSolver::UpdateVectorSolution, this, std::ref(Keff), std::ref(Phi),
				std::ref(B_r), ie);

			// Compute the frequency-dependent stiffness matrix
			std::thread thread1(&JacobiDavidsonNLEigenSolver::getFreqDependentStiffMtx, this, std::ref(K0), std::ref(MM),
				std::ref(Kn), std::ref(Omega(ie)));

			// Compute the frequency-dependent mass matrix
			std::thread thread2(&JacobiDavidsonNLEigenSolver::getFreqDependentMassMtx, this, std::ref(MM), std::ref(Mn), std::ref(Omega(ie)));
			thread0.join();
			thread1.join();
			thread2.join();

			PtMP = Phi.col(ie).transpose() * Mn * Phi.col(ie);
			PtKP = Phi.col(ie).transpose() * Kn * Phi.col(ie);
			theta = PtKP / PtMP;

			//LOG_ASSERT(!(PtMP < 0), "Error: Negative mass matrix!!!");

			// Normalize the improved Vector
			Phi.col(ie) = (1.0 / sqrt(PtMP)) * Phi.col(ie);

			// Evaluate the convergence
			conv = (theta - Omega(ie)) / theta;
			LOG_ASSERT(!isnan(conv), "Error: Not-a-number in the computed eigenvalues!");

			LOG_INFO("iter: {0}    rel.error: {1}", iterK, abs((double)conv));

			//Update the new eigenvalue
			Omega(ie) = theta;

			// Check the max. number of iterations
			iterK++;

			if (iterK > m_MaxIter)
			{
				LOG_ASSERT(false, "Error: It has reached the max. number of iterations!!");
				return false;
			}

		}
	}

	//Print results
	printResults(Omega, Phi);

	LOG_INFO("NLEigenvalue routine has run sucessfully!");
	return true;
}

void JacobiDavidsonNLEigenSolver::readFileAndGetStiffMassMatrices(SparseMatrix& K0, std::vector<SparseMatrix>& MM, Vector& Omega)
{
	//Open file to read
	std::fstream fid;
	fid.open(m_FilePath, std::ios::in);

	//Check error
	if (fid.fail())
	{
		LOG_ASSERT(fid.fail(), "ERROR: Error in opening the file!");
	}



	if (fid.is_open())
	{
		std::string line;
		std::getline(fid, line);
		// Read #dof, #mass matrices, #eigenvalues
		int hasEigenTrial;
		fid >> m_Dimensions >> m_NumberOfMassMtx >> m_NumberOfEigenValues >> m_TOL >> hasEigenTrial;

		if (hasEigenTrial == 1)
		{
			m_hasInitialTrial = true;
		}

		// Set the matrices
		K0 = SparseMatrix(m_Dimensions, m_Dimensions);
		Omega = Vector(m_NumberOfEigenValues);
		Omega.setZero(); 
		MM.reserve(m_NumberOfMassMtx);

		// Read the stiffness matrix K0
		for (int ii = 0; ii < m_Dimensions; ii++)
		{
			for (int jj = 0; jj < m_Dimensions; jj++)
			{
				double valueIJ;
				fid >> valueIJ;
				// It stores only the nonzero elements into the K0 matrix
				if (!(valueIJ == 0.0))
					K0.insert(ii, jj) = valueIJ;
			    
			}
		}

		// Itsuppresses the remaining empty spaceand transforms the matrix into
		// a compressed column storage
		K0.makeCompressed();

		//Read mass matrices
		for (int im = 0; im < m_NumberOfMassMtx; im++)
		{
			SparseMatrix Mtemp(m_Dimensions, m_Dimensions);
			for (int ii = 0; ii < m_Dimensions; ii++)
			{
				for (int jj = 0; jj < m_Dimensions; jj++)
				{
					// Get value
					double valueIJ;
					fid >> valueIJ;
					// It stores only the nonzero elements into the K0 matrix
					if (!(valueIJ == 0.0))
						Mtemp.insert(ii, jj) = valueIJ;
					
				}
			}

			Mtemp.makeCompressed();
			MM.emplace_back(-Mtemp);
		}

		// Check 
		if (m_hasInitialTrial)
		{
			LOG_INFO("A first attemptive of the eigenvalues has been provided!");
			for (int ie = 0; ie < m_NumberOfEigenValues; ie++)
			{
				fid >> Omega(ie);
			}
		}


	}

	//Close file
	fid.close();
}

void JacobiDavidsonNLEigenSolver::printResults(Vector& Omega, DenseMatrix& Phi) const
{
	// Save the eigenproblem results
	LOG_INFO("Save the eigenvalues and Vector!");

	// Get the directory path
	std::string directory, resultFile1, resultFile2;
	const size_t last_slash_idx = m_FilePath.rfind('\\');
	if (std::string::npos != last_slash_idx)
	{
		directory = m_FilePath.substr(0, last_slash_idx);
	}

	resultFile1 = directory + "\\Phi.dat";
	resultFile2 = directory + "\\Omega.dat";

	std::ofstream out1, out2;
	out1.open(resultFile1);
	out2.open(resultFile2);

	if (!out1 || !out2)
	{
		LOG_ASSERT(false, "ERROR: Error in opening the file!");
	}

	//Save Phi
	//out1 << m_Dimensions <<  " " << m_NumberOfEigenValues << std::endl;
	out1 << std::setprecision(16) << std::scientific << Phi;
	out1.close();

	//Save Omega
	//out2 << m_NumberOfEigenValues << std::endl;
	out2 << std::setprecision(16) << std::scientific << Omega;
	out2.close();
}

void JacobiDavidsonNLEigenSolver::UpdateVectorSolution(DenseMatrix& Keff, DenseMatrix& Phi, DenseMatrix& B_r, int index)
{
	Vector rk, dUk(m_Dimensions);

	// Evaluate the residual error
	rk = -Keff * Phi.col(index);
	dUk.setOnes();

	// Project the effective stiffness matrix
	projectEffectiveStiffMatrix(Keff, B_r, index);

	// solve dUk
	iterativeLinearSolver(Keff, rk, dUk);

	//Project dUk
	for (int is = 0; is < index + 1; is++)
	{
		dUk += -B_r.col(is) * (B_r.col(is).transpose() * dUk);
	}

	// Update solution
	Phi.col(index) += dUk;
}

void JacobiDavidsonNLEigenSolver::getFreqDependentStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, SparseMatrix& Kn, data_type omega)
{
	//Initialize
	Kn = K0;

	for (int jj = 1; jj < m_NumberOfMassMtx; jj++)
	{

		Kn += (jj)*pow(omega, jj + 1.0) * MM[jj];
	}
}

void JacobiDavidsonNLEigenSolver::getFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mn, data_type omega)
{
	//Initialize
	Mn = MM[0];

	for (int jj = 1; jj < m_NumberOfMassMtx; jj++)
	{

		Mn += (jj + 1.0) * pow(omega, jj) * MM[jj];
	}
}

void JacobiDavidsonNLEigenSolver::getGeneralizedFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mlrls, data_type lr, data_type ls)
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

void JacobiDavidsonNLEigenSolver::getEffectiveStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, DenseMatrix& Keff, data_type omega)
{
	// Initialize
	Keff = K0;

	for (int jj = 0; jj < m_NumberOfMassMtx; jj++)
	{
		Keff -= pow(omega, jj + 1.0) * MM[jj];
	}
}

void JacobiDavidsonNLEigenSolver::projectEffectiveStiffMatrix(DenseMatrix& Keff, DenseMatrix& B_s, int indexEig)
{
	// Project the effective stiffness matrix onto the subspace
	// orthogonal to all preceding Vectors and add the orthogonal
	// projector to make it nonsingular
	for (int ii = 0; ii < indexEig + 1; ii++)
	{
		Keff += (B_s.col(ii) - Keff * B_s.col(ii)) * (B_s.col(ii).transpose());
	
	}

}

bool JacobiDavidsonNLEigenSolver::iterativeLinearSolver(const DenseMatrix& A, const Vector& b, Vector& x)
{
	// Implementation of the non-preconditioned BICGSTAB algorithm
	//Vector r(m_Dimensions), r_hat(m_Dimensions);
	//auto bmod = A.transpose() * b;
	//auto Amod = A.transpose() * A;
	//r = bmod - Amod*x;
	//r_hat = r;
	//Vector v(m_Dimensions), p(m_Dimensions), t(m_Dimensions), s(m_Dimensions), h(m_Dimensions);
	//v.setZero(); p.setZero(); t.setZero(); s.setZero(); h.setZero();
	//double rho = 1.0, alpha = 1.0, omega = 1.0, beta, tol = 1e-10;
	//double error, bnorm2 = bmod.norm();
	//int icount, max_iter = 2000;
	bool status = true;

	//for (int i = 0; i < max_iter; i++)
	//{
	//	beta = 1.0 / rho * (alpha / omega);
	//	rho = r_hat.transpose() * r;
	//	beta *= beta;

	//	p = r + beta * (p - omega * v);
	//	v = Amod * p;
	//	alpha = rho / (r_hat.transpose()*v);
	//	h = x + alpha * p;

	//	s = r - alpha * v;
	//	t = Amod * s;
	//	omega = (t.transpose() * s);
	//	omega *= 1.0 / (t.transpose() * t);
	//	
	//	// Update solution
	//	x = h + omega * s;

	//	// Check solution
	//	error = (bmod - Amod * x).norm();
	//	error *= 1.0 / bnorm2;
	//	if (error < tol)
	//	{
	//		icount = i;
	//		break;
	//	}

	//	r = s - omega * t;
	//}

	//if (error > tol)
	//{
	//	status = false;
	//	icount = max_iter;
	//	
	//}
	// Direct solve
	//auto bmod = A.transpose() * b;
	//auto Amod = A.transpose() * A;
	//x = Amod.fullPivLu().solve(bmod);
	//
	// 
	//data_type bnorm2 = bmod.norm();
	//data_type error = (bmod - Amod * x).norm() / (bnorm2);

	x = A.fullPivLu().solve(b);
	//data_type bnorm2 = b.norm();
	//data_type error = (b - A * x).norm() / (bnorm2);

	// Check the L2 norm
	//auto res = (A * x - b);
	//data_type relative_error1 = (A * x - b).norm();
	//LOG_WARN("Relative error(L2 norm) BDCSVD = {0}", relative_error1);
	//LOG_WARN("Error(L2 norm) = {0}", error);

	//LOG_INFO("#Iterations: {0},    Estimated error: {1}", linsolver.iterations(), linsolver.error());



	/*if (linsolver.error() > 1.0e-14)
	{
		status = false;
		LOG_ASSERT(status, "It has reach the max number of iterations!");
	}
		*/
	return status;
}


