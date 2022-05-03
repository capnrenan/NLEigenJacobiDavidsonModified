#include "nlpch.h"
#include "NLEigenJacobiDavidson.h"

#ifdef ENABLE_BLAZE

// Method implementations using Eigen
NLEigenJacobiDavidson::NLEigenJacobiDavidson(const std::string& filepath)
	: m_Dimensions(0), m_NumberOfMassMtx(1), m_NumberOfEigenValues(0),
	m_MaxIter(20), m_TOL(1e-14), m_FilePath(filepath)
{
	// Initialize the log system
	Log::Init();
	blaze::setNumThreads(8);
}

NLEigenJacobiDavidson::~NLEigenJacobiDavidson()
{

}

bool NLEigenJacobiDavidson::execute()
{
	//Reading data
	LOG_INFO("Reading filedata...\n");
	blaze::DynamicMatrix<double, blaze::rowMajor> K0;
	std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>> MM;
	readFileAndGetStiffMassMatrices(K0, MM);

	LOG_INFO("Initializing the matrices...\n");
	//Initialize the matrices and set them as zero
	blaze::DynamicVector<double, blaze::columnVector> Omega(m_NumberOfEigenValues);
	blaze::DynamicVector<double, blaze::columnVector> rk(m_Dimensions), dUk(m_Dimensions);
	blaze::DynamicMatrix<double, blaze::rowMajor> Phi(m_Dimensions, m_NumberOfEigenValues);
	blaze::DynamicMatrix<double, blaze::rowMajor> B_r(m_Dimensions, m_NumberOfEigenValues);
	blaze::DynamicMatrix<double, blaze::rowMajor> Keff(m_Dimensions, m_Dimensions);
	blaze::DynamicMatrix<double, blaze::rowMajor> Kn(m_Dimensions, m_Dimensions);
	blaze::DynamicMatrix<double, blaze::rowMajor> Mn(m_Dimensions, m_Dimensions);
	blaze::DynamicMatrix<double, blaze::rowMajor> Mlrls(m_Dimensions, m_Dimensions);

	LOG_INFO("Solving with TOL = {0}\n", m_TOL);

	//Set as zero
	Omega = 0.0;   rk = 0.0; dUk = 0.0; B_r = 0.0; Phi = 1.0;
	Keff = 0.0; Kn = 0.0; Mn = 0.0; Mlrls = 0.0;

	//Auxiliary variables
	double conv, normBr, PtMP, PtKP, theta;
	int iterK;

	LOG_INFO("Processing...");
	// Loop in each eigenvalue
	for (int ie = 0; ie < m_NumberOfEigenValues; ie++)
	{
		// Set the convergence parameters
		conv = 1.0;
		iterK = 0;

		LOG_INFO("---------------------------------\nEigenvalue #{0}:", ie);

		if (ie > 0)
		{
			Omega[ie] = Omega[ie - 1];
		}

		while (abs(conv) > m_TOL)
		{
			// Orthogonalized phi_r with respect to phi_s
			for (int is = 0; is < ie + 1; is++)
			{
				// Get the generalized freq-dependent mass matrix
				getGeneralizedFreqDependentMassMtx(MM, Mlrls, Omega[ie], Omega[is]);
				blaze::column(B_r, ie) = Mlrls * blaze::column(Phi, ie);
				//B_r.col(ie) = Mlrls * Phi.col(is);


				if (is > 0)
				{
					for (int el = 0; el < is; el++)
					{
						//Project b_s = b_s - b_el*(b_el.t()*b_s)
						//B_r.col(is) += -B_r.col(el) * (B_r.col(el).transpose() * B_r.col(is));
						blaze::column(B_r,is) += -blaze::column(B_r, el) * (blaze::trans(blaze::column(B_r, el))*blaze::column(B_r, is));
					}
				}

				// Normalize
				//normBr = sqrt(B_r.col(is).transpose() * B_r.col(is));
				//B_r.col(is) = 1.0 / (normBr)*B_r.col(is);
				normBr = std::sqrt(blaze::trans(blaze::column(B_r,is)) * blaze::column(B_r, is));
				blaze::column(B_r, is) *= (1.0 / normBr);
		
			}

			// Orthogonalize phi_e with respect to the preceding eigenvector phi
			if (ie > 0)
			{
				for (int is = 0; is < ie; is++)
				{
					//Phi.col(ie) += -B_r.col(is) * (B_r.col(is).transpose() * Phi.col(ie));
					blaze::column(Phi, ie) -= blaze::column(B_r, is) * (blaze::column(B_r, is), blaze::column(Phi, ie));
				}
			}

			// Evaluate the effective stiffness matrix
			getEffectiveStiffMtx(K0, MM, Keff, Omega[ie]);

			// Evaluate the residual error
			rk = -Keff * blaze::column(Phi,ie);


			// Project the effective stiffness matrix
			projectEffectiveStiffMatrix(Keff, B_r, ie);

			// solve dUk
			//std::cout << rk << std::endl;
			iterativeLinearSolver(Keff, rk, dUk);

			//Project dUk
			for (int is = 0; is < ie + 1; is++)
			{
				//dUk += -B_r.col(is) * (B_r.col(is).transpose() * dUk);
				dUk += -blaze::column(B_r,is) * (blaze::column(B_r, is), dUk);
			}

			// Update solution
			//Phi.col(ie) += dUk;
			blaze::column(Phi,ie) += dUk;

			// Evaluate the Rayleigh quotient
			getFreqDependentStiffMtx(K0, MM, Kn, Omega[ie]);
			getFreqDependentMassMtx(MM, Mn, Omega[ie]);
			//PtMP = Phi.col(ie).transpose() * Mn * Phi.col(ie);
			//PtKP = Phi.col(ie).transpose() * Kn * Phi.col(ie);
			PtMP = (blaze::column(Phi, ie), Mn * blaze::column(Phi, ie));
			PtKP = (blaze::column(Phi, ie), Kn * blaze::column(Phi, ie));
			theta = PtKP / PtMP;

			LOG_ASSERT(!(PtMP < 0), "Error: Negative mass matrix!!!");

			// Normalize the improved eigenvector
			//Phi.col(ie) = (1.0 / sqrt(PtMP)) * Phi.col(ie);
			blaze::column(Phi, ie) = (1.0 / sqrt(PtMP)) * blaze::column(Phi, ie);

			// Evaluate the convergence
			conv = abs(theta - Omega[ie]) / theta;

			LOG_INFO("iter: {0}    rel.error: {1}", iterK, conv);

			//Update the new eigenvalue
			Omega[ie] = theta;

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

	return true;
}

void NLEigenJacobiDavidson::readFileAndGetStiffMassMatrices(blaze::DynamicMatrix<double, blaze::rowMajor>& K0, std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM)
{
	//Open file to read
	std::fstream fid;
	fid.open(m_FilePath, std::ios::in);

	//Check error
	if (fid.fail())
	{
		LOG_ASSERT(false, "ERROR: Error in opening the file!");
	}

	if (fid.is_open())
	{
		std::string line;
		std::getline(fid, line);
		// Read #dof, #mass matrices, #eigenvalues
		fid >> m_Dimensions >> m_NumberOfMassMtx >> m_NumberOfEigenValues >> m_TOL;

		// Set the matrices
		K0 = blaze::DynamicMatrix<double, blaze::rowMajor>(m_Dimensions, m_Dimensions);
		K0 = 0.0;
		blaze::DynamicMatrix<double, blaze::rowMajor> Mtemp(m_Dimensions, m_Dimensions);
		Mtemp = 0.0;
		MM.reserve(m_NumberOfMassMtx);

		// Read the stiffness matrix K0
		for (int ii = 0; ii < m_Dimensions; ii++)
		{
			for (int jj = 0; jj < m_Dimensions; jj++)
			{
				fid >> K0(ii, jj);
			}
		}

		//Read mass matrices
		for (int im = 0; im < m_NumberOfMassMtx; im++)
		{
			for (int ii = 0; ii < m_Dimensions; ii++)
			{
				for (int jj = 0; jj < m_Dimensions; jj++)
				{
					fid >> Mtemp(ii, jj);
				}
			}

			MM.emplace_back(-Mtemp);
		}
	}

	//Close file
	fid.close();
}

void NLEigenJacobiDavidson::printResults(blaze::DynamicVector<double, blaze::columnVector>& Omega, blaze::DynamicMatrix<double, blaze::rowMajor>& Phi) const
{
	// Save the eigenproblem results
	LOG_INFO("Save the eigenvalues and eigenvector!");

	// Get the directory path
	std::string directory, resultFile1, resultFile2;
	const size_t last_slash_idx = m_FilePath.rfind('\/');
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
	//out1 << std::setprecision(16) << std::scientific << Phi;
	for (int ii = 0; ii < Phi.columns(); ii++)
	{
		for (int jj = 0; jj < Phi.rows(); jj++)
		{
			out1 << std::setprecision(16) << std::scientific << Phi(ii,jj) << " ";
		}
		out1 << std::endl;
	}
	out1.close();

	//Save Omega
	//out2 << m_NumberOfEigenValues << std::endl;
	for (int ii = 0; ii < Omega.size(); ii++)
	{
		out2 << std::setprecision(16) << std::scientific << Omega[ii] << std::endl;
	}
	//out2 << std::setprecision(16) << std::scientific << Omega;
	out2.close();
}

void NLEigenJacobiDavidson::getFreqDependentStiffMtx(const blaze::DynamicMatrix<double, blaze::rowMajor>& K0, const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Kn, double omega)
{
	//Initialize
	Kn = K0;

	for (int jj = 1; jj < m_NumberOfMassMtx; jj++)
	{
		Kn += (jj)*pow(omega, jj + 1.0) * MM[jj];
	}
}

void NLEigenJacobiDavidson::getFreqDependentMassMtx(const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Mn, double omega)
{
	//Initialize
	Mn = MM[0];

	for (int jj = 1; jj < m_NumberOfMassMtx; jj++)
	{

		Mn += (jj + 1.0) * pow(omega, jj) * MM[jj];
	}
}

void NLEigenJacobiDavidson::getGeneralizedFreqDependentMassMtx(const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Mlrls, double lr, double ls)
{
	//Initialize
	Mlrls = 0.0;
	for (int jj = 0; jj < m_NumberOfMassMtx; jj++)
	{
		for (int kk = 0; kk < jj + 1; kk++)
		{
			Mlrls += pow(lr, kk) * pow(ls, jj - kk) * MM[jj];
		}
	}
}

void NLEigenJacobiDavidson::getEffectiveStiffMtx(const blaze::DynamicMatrix<double, blaze::rowMajor>& K0, const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Keff, double omega)
{
	// Initialize
	Keff = K0;

	for (int jj = 0; jj < m_NumberOfMassMtx; jj++)
	{
		Keff -= pow(omega, jj + 1.0) * MM[jj];
	}
}

void NLEigenJacobiDavidson::projectEffectiveStiffMatrix(blaze::DynamicMatrix<double, blaze::rowMajor>& Keff, blaze::DynamicMatrix<double, blaze::rowMajor>& B_s, int indexEig)
{
	// Project the effective stiffness matrix onto the subspace
	// orthogonal to all preceding eigenvectors and add the orthogonal
	// projector to make it nonsingular
	for (int ii = 0; ii < indexEig + 1; ii++)
	{
		Keff += (blaze::column(B_s,ii) - Keff * blaze::column(B_s, ii)) * blaze::trans(blaze::column(B_s, ii));
	}
}

bool NLEigenJacobiDavidson::iterativeLinearSolver(const blaze::DynamicMatrix<double, blaze::rowMajor>& A, const blaze::DynamicVector<double, blaze::columnVector>& b, blaze::DynamicVector<double, blaze::columnVector>& x)
{
	//auto temp  = blaze::solve(A, b);
	//x = temp;
	
	double error;
	bool status = true;
// #pragma omp parallel
	{
		blaze::solve(A, x, b);   // Computing the solution x
		auto temp = A * x - b;
		error = (temp, temp);
	}

	if (error > m_TOL)
		status = false;
	
	return status;
}


#else

// Method implementations using Eigen
NLEigenJacobiDavidson::NLEigenJacobiDavidson(const std::string& filepath)
	: m_Dimensions(0), m_NumberOfMassMtx(1), m_NumberOfEigenValues(0),
	m_MaxIter(20), m_TOL(1e-14), m_FilePath(filepath)
{
	// Initialize the log system
	Log::Init();
	Eigen::initParallel();
}

NLEigenJacobiDavidson::~NLEigenJacobiDavidson()
{

}

bool NLEigenJacobiDavidson::execute()
{
	//Reading data
	LOG_INFO("Reading filedata...\n");
	Eigen::MatrixXd K0;
	std::vector<Eigen::MatrixXd> MM;
	readFileAndGetStiffMassMatrices(K0, MM);

	LOG_INFO("Initializing the matrices...\n");
	//Initialize the matrices and set them as zero
	Eigen::VectorXd Omega(m_NumberOfEigenValues);
	Eigen::VectorXd rk, dUk;
	Eigen::MatrixXd Phi(m_Dimensions, m_NumberOfEigenValues);
	Eigen::MatrixXd B_r(m_Dimensions, m_NumberOfEigenValues);
	Eigen::MatrixXd Keff(m_Dimensions, m_Dimensions);
	Eigen::MatrixXd Kn(m_Dimensions, m_Dimensions);
	Eigen::MatrixXd Mn(m_Dimensions, m_Dimensions);
	Eigen::MatrixXd Mlrls(m_Dimensions, m_Dimensions);

	LOG_INFO("Solving with TOL = {0}\n", m_TOL);

	//Set as zero
	Omega.setZero();   B_r.setZero(); Phi.setOnes();
	Keff.setZero(); Kn.setZero(); Mn.setZero(); Mlrls.setZero();

	//Auxiliary variables
	double conv, normBr, PtMP, PtKP, theta;
	int iterK;

	LOG_INFO("Processing...");
	// Loop in each eigenvalue
	for (int ie = 0; ie < m_NumberOfEigenValues; ie++)
	{
		// Set the convergence parameters
		conv = 1.0;
		iterK = 0;

		LOG_INFO("---------------------------------\nEigenvalue #{0}:", ie);

		if (ie > 0)
		{
			Omega(ie) = Omega(ie - 1);
		}

		while (abs(conv) > m_TOL)
		{
			// Orthogonalized phi_r with respect to phi_s
			for (int is = 0; is < ie + 1; is++)
			{
				// Get the generalized freq-dependent mass matrix
				getGeneralizedFreqDependentMassMtx(MM, Mlrls, Omega(ie), Omega(is));
				B_r.col(ie) = Mlrls * Phi.col(is);

				if (is > 0)
				{
					for (int el = 0; el < is; el++)
					{
						//Project b_s = b_s - b_el*(b_el.t()*b_s)
						B_r.col(is) += -B_r.col(el) * (B_r.col(el).transpose() * B_r.col(is));
					}
				}

				// Normalize
				normBr = sqrt(B_r.col(is).transpose() * B_r.col(is));
				B_r.col(is) = 1.0 / (normBr)*B_r.col(is);
			}

			// Orthogonalize phi_e with respect to the preceding eigenvector phi
			if (ie > 0)
			{
				for (int is = 0; is < ie; is++)
				{
					Phi.col(ie) += -B_r.col(is) * (B_r.col(is).transpose() * Phi.col(ie));
				}
			}

			// Evaluate the effective stiffness matrix
			getEffectiveStiffMtx(K0, MM, Keff, Omega(ie));

			// Evaluate the residual error
			rk = -Keff * Phi.col(ie);

			// Project the effective stiffness matrix
			projectEffectiveStiffMatrix(Keff, B_r, ie);

			// solve dUk
			iterativeLinearSolver(Keff, rk, dUk);

			//Project dUk
			for (int is = 0; is < ie + 1; is++)
			{
				dUk += -B_r.col(is) * (B_r.col(is).transpose() * dUk);
			}

			// Update solution
			Phi.col(ie) += dUk;

			// Evaluate the Rayleigh quotient
			getFreqDependentStiffMtx(K0, MM, Kn, Omega(ie));
			getFreqDependentMassMtx(MM, Mn, Omega(ie));
			PtMP = Phi.col(ie).transpose() * Mn * Phi.col(ie);
			PtKP = Phi.col(ie).transpose() * Kn * Phi.col(ie);
			theta = PtKP / PtMP;

			LOG_ASSERT(!(PtMP < 0), "Error: Negative mass matrix!!!");

			// Normalize the improved eigenvector
			Phi.col(ie) = (1.0 / sqrt(PtMP)) * Phi.col(ie);

			// Evaluate the convergence
			conv = abs(theta - Omega(ie)) / theta;

			LOG_INFO("iter: {0}    rel.error: {1}", iterK, conv);

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

	return true;
}


void NLEigenJacobiDavidson::readFileAndGetStiffMassMatrices(Eigen::MatrixXd& K0, std::vector<Eigen::MatrixXd>& MM)
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
		fid >> m_Dimensions >> m_NumberOfMassMtx >> m_NumberOfEigenValues >> m_TOL;

		// Set the matrices
		K0 = Eigen::MatrixXd(m_Dimensions, m_Dimensions);
		Eigen::MatrixXd Mtemp(m_Dimensions, m_Dimensions);
		K0.setZero(); Mtemp.setZero();
		MM.reserve(m_NumberOfMassMtx);

		// Read the stiffness matrix K0
		for (int ii = 0; ii < m_Dimensions; ii++)
		{
			for (int jj = 0; jj < m_Dimensions; jj++)
			{
				fid >> K0(ii, jj);
			}
		}

		//Read mass matrices
		for (int im = 0; im < m_NumberOfMassMtx; im++)
		{
			for (int ii = 0; ii < m_Dimensions; ii++)
			{
				for (int jj = 0; jj < m_Dimensions; jj++)
				{
					fid >> Mtemp(ii, jj);
				}
			}

			MM.emplace_back(-Mtemp);
		}
	}

	//Close file
	fid.close();
}

void NLEigenJacobiDavidson::printResults(Eigen::VectorXd& Omega, Eigen::MatrixXd& Phi) const
{
	// Save the eigenproblem results
	LOG_INFO("Save the eigenvalues and eigenvector!");

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
	for (int ii = 0; ii < indexEig + 1; ii++)
	{
		Keff += (B_s.col(ii) - Keff * B_s.col(ii)) * (B_s.col(ii).transpose());
	}
}

bool NLEigenJacobiDavidson::iterativeLinearSolver(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{
	// Set the iterative linear solver (Conjugate Gradients)
	// The number of max. of iterations
	//Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower> linsolver;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> linsolver;
	//Eigen::BiCGSTAB<Eigen::MatrixXd> linsolver;
	//linsolver.setTolerance(1.0e-14);
	//linsolver.setMaxIterations(4 * m_Dimensions);

	//Solve
	//x = linsolver.compute(A.transpose()*A).solve(A.transpose()*b);

	// Direct solve
	x = A.fullPivLu().solve(b);
	//Eigen::LLT<Eigen::MatrixXd> llt;

	//llt.compute(A);
	//x = llt.solve(b);
	//x = (A.adjoint() * A).llt().solve((A.adjoint() * b));
	//x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);


	// Check the L2 norm
	double relative_error1 = (A * x - b).norm() / b.norm();

	//LOG_WARN("Relative error(L2 norm) BDCSVD = {0}", relative_error1);
	//LOG_WARN("Relative error(L2 norm) fullPivLu = {0}", relative_error1);

	//LOG_INFO("#Iterations: {0},    Estimated error: {1}", linsolver.iterations(), linsolver.error());

	bool status = true;

	/*if (linsolver.error() > 1.0e-14)
	{
		status = false;
		LOG_ASSERT(status, "It has reach the max number of iterations!");
	}
		*/
	return status;
}


#endif
