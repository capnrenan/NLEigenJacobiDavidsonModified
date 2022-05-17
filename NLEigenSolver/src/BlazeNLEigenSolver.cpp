#include "nlpch.h"
#include "BlazeNLEigenSolver.h"

BlazeNLEigenSolver::BlazeNLEigenSolver(const std::string& filepath)
	: m_Dimensions(0), m_NumberOfMassMtx(1), m_NumberOfEigenValues(0),
	m_MaxIter(20), m_TOL(1e-14), m_FilePath(filepath)
{
	blaze::setNumThreads(10);
}

BlazeNLEigenSolver::~BlazeNLEigenSolver()
{
}

bool BlazeNLEigenSolver::execute()
{
	PROFILE_FUNCTION();
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
						blaze::column(B_r, is) += -blaze::column(B_r, el) * (blaze::trans(blaze::column(B_r, el)) * blaze::column(B_r, is));
					}
				}

				// Normalize
				//normBr = sqrt(B_r.col(is).transpose() * B_r.col(is));
				//B_r.col(is) = 1.0 / (normBr)*B_r.col(is);
				normBr = std::sqrt(blaze::trans(blaze::column(B_r, is)) * blaze::column(B_r, is));
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
			rk = -Keff * blaze::column(Phi, ie);


			// Project the effective stiffness matrix
			projectEffectiveStiffMatrix(Keff, B_r, ie);

			// solve dUk
			//std::cout << rk << std::endl;
			iterativeLinearSolver(Keff, rk, dUk);

			//Project dUk
			for (int is = 0; is < ie + 1; is++)
			{
				//dUk += -B_r.col(is) * (B_r.col(is).transpose() * dUk);
				dUk += -blaze::column(B_r, is) * (blaze::column(B_r, is), dUk);
			}

			// Update solution
			//Phi.col(ie) += dUk;
			blaze::column(Phi, ie) += dUk;

			// Evaluate the Rayleigh quotient
			//getFreqDependentStiffMtx(K0, MM, Kn, Omega[ie]);
			//getFreqDependentMassMtx(MM, Mn, Omega[ie]);

			std::thread thread(&BlazeNLEigenSolver::getFreqDependentStiffMtx, this, std::ref(K0), std::ref(MM),
				std::ref(Kn), std::ref(Omega[ie]));
			std::thread thread2(&BlazeNLEigenSolver::getFreqDependentMassMtx, this, std::ref(MM), std::ref(Mn), std::ref(Omega[ie]));
			thread.join();
			thread2.join();

			//PtMP = Phi.col(ie).transpose() * Mn * Phi.col(ie);
			//PtKP = Phi.col(ie).transpose() * Kn * Phi.col(ie);
			PtMP = (blaze::column(Phi, ie), Mn * blaze::column(Phi, ie));
			PtKP = (blaze::column(Phi, ie), Kn * blaze::column(Phi, ie));
			/*#pragma omp sections
			{
				#pragma omp section
				PtMP = (blaze::column(Phi, ie), Mn * blaze::column(Phi, ie));

				#pragma omp section
				PtKP = (blaze::column(Phi, ie), Kn * blaze::column(Phi, ie));
			}*/

			//PtMP = (blaze::column(Phi, ie), Mn * blaze::column(Phi, ie));
			
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
	//printResults(Omega, Phi);

	return true;
}

void BlazeNLEigenSolver::readFileAndGetStiffMassMatrices(blaze::DynamicMatrix<double, blaze::rowMajor>& K0, std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM)
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

void BlazeNLEigenSolver::printResults(blaze::DynamicVector<double, blaze::columnVector>& Omega, blaze::DynamicMatrix<double, blaze::rowMajor>& Phi) const
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
			out1 << std::setprecision(16) << std::scientific << Phi(ii, jj) << " ";
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

void BlazeNLEigenSolver::getFreqDependentStiffMtx(const blaze::DynamicMatrix<double, blaze::rowMajor>& K0, const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Kn, double omega)
{
	//Initialize
	Kn = K0;

	for (int jj = 1; jj < m_NumberOfMassMtx; jj++)
	{
		Kn += (jj)*pow(omega, jj + 1.0) * MM[jj];
	}
}

void BlazeNLEigenSolver::getFreqDependentMassMtx(const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Mn, double omega)
{
	//Initialize
	Mn = MM[0];

	for (int jj = 1; jj < m_NumberOfMassMtx; jj++)
	{

		Mn += (jj + 1.0) * pow(omega, jj) * MM[jj];
	}
}

void BlazeNLEigenSolver::getGeneralizedFreqDependentMassMtx(const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Mlrls, double lr, double ls)
{
	PROFILE_FUNCTION();
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

void BlazeNLEigenSolver::getEffectiveStiffMtx(const blaze::DynamicMatrix<double, blaze::rowMajor>& K0, const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Keff, double omega)
{
	PROFILE_FUNCTION();
	// Initialize
	Keff = K0;

	for (int jj = 0; jj < m_NumberOfMassMtx; jj++)
	{
		Keff -= pow(omega, jj + 1.0) * MM[jj];
	}
}

void BlazeNLEigenSolver::projectEffectiveStiffMatrix(blaze::DynamicMatrix<double, blaze::rowMajor>& Keff, blaze::DynamicMatrix<double, blaze::rowMajor>& B_s, int indexEig)
{
	// Project the effective stiffness matrix onto the subspace
	// orthogonal to all preceding eigenvectors and add the orthogonal
	// projector to make it nonsingular
	for (int ii = 0; ii < indexEig + 1; ii++)
	{
		Keff += (blaze::column(B_s, ii) - Keff * blaze::column(B_s, ii)) * blaze::trans(blaze::column(B_s, ii));
	}
}

bool BlazeNLEigenSolver::iterativeLinearSolver(const blaze::DynamicMatrix<double, blaze::rowMajor>& A, const blaze::DynamicVector<double, blaze::columnVector>& b, blaze::DynamicVector<double, blaze::columnVector>& x)
{
	//auto temp  = blaze::solve(A, b);
	//x = temp;
	PROFILE_FUNCTION();
	double error;
	bool status = true;
	{
		blaze::solve(A, x, b);   // Computing the solution x
		auto temp = A * x - b;
		error = (temp, temp);
	}

	if (error > m_TOL)
		status = false;

	return status;
}
