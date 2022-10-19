#include "nlpch.h"
#include "inverseFreeKrylovNLEigenSolver.h"

// Lapack wrappers
#include "f2c.h"
#include "clapack.h"

inverseFreeKrylovNLEigenSolver::inverseFreeKrylovNLEigenSolver(const std::string& filepath)
    : m_FilePath(filepath)
{


}

inverseFreeKrylovNLEigenSolver::~inverseFreeKrylovNLEigenSolver()
{

}

bool inverseFreeKrylovNLEigenSolver::execute()
{
	///Reading data
	std::cout << "Reading filedata: " << m_FilePath << std::endl << std::endl;
	//LOG_INFO("Reading filedata: {0}\n", m_FilePath);
	SparseMatrix K0;
	std::vector<SparseMatrix> MM_0;
	std::vector<DenseMatrix>MM_Update;
	readFileAndGetStiffMassMatrices(K0, MM_0);
	DenseMatrix Kt = K0;

	// Print info about the problem
	printInfo();

	// Set the MM_Update as MM_0
	MM_Update.reserve(m_InputData.NumberOfMassMtx);
	for (int im = 0; im < m_InputData.NumberOfMassMtx; im++)
	{
		MM_Update.emplace_back(MM_0[im]);
	}
	
	// Initialize the eigenvalues and eigenvector to compute
	Vector Omega(m_InputData.NumberOfEigenvalues);
	DenseMatrix Phi(m_InputData.Dimensions, m_InputData.NumberOfEigenvalues);
	DenseMatrix Mlrls(m_InputData.Dimensions, m_InputData.Dimensions);
	Omega.setZero(); Phi.setOnes(); Mlrls.setZero();

	// initialize status
	resultStatus resultStatus;

	// set temp values
	data_type omega_temp;
	Vector phi_temp;

	// Loop to compute the k smallest eigenpairs of (K(omega),M(omega))
	for (int k = 0; k < m_InputData.NumberOfEigenvalues; k++)
	{
		LOG_INFO("---------------------------------\nEigenvalue #{0}:", k);

		if (k > 0)
		{
			// Get the previous compute eigenpair
			omega_temp = Omega(k - 1);
			phi_temp = Phi.col(k - 1);

			// DEFLATION PROCEDURE
			generalizedDeflationProcedure(Kt, MM_Update, Mlrls, omega_temp, phi_temp);

			// Compute the smallest eigenpair of the deflated pencil
			computeSmallestEigenpairKrylov(Kt, MM_Update, MM_0, omega_temp, phi_temp, 
											m_InputData.NumberOfKrylovBasis, resultStatus);


			// Set the computed eigenpairs
			Omega(k) = omega_temp;
			Phi.col(k) = phi_temp;

		}
		else
		{
			omega_temp = Omega(k);
			phi_temp = Phi.col(k);
			computeSmallestEigenpairKrylov(K0, MM_Update, MM_0, omega_temp, phi_temp,
											m_InputData.NumberOfKrylovBasis, resultStatus);

			Omega(k) = omega_temp;
			Phi.col(k) = phi_temp;
		}

		LOG_INFO("#Iter: {0}, conv = {1}", resultStatus.numberOfIterations, resultStatus.Convergence);

		if (!resultStatus.status)
		{
			LOG_WARN("Error: It has not converged!");
			return false;
		}
	}
	
	LOG_INFO("--------------------------------------------------");
	LOG_INFO("Saving eigenpair solutions in file!");
	printResults(Omega, Phi);
	LOG_INFO("Done!");


	return true;
}

bool inverseFreeKrylovNLEigenSolver::findEigenvaluesFromInitialGuess()
{
	return false;
}

void inverseFreeKrylovNLEigenSolver::readFileAndGetStiffMassMatrices(SparseMatrix& K0, std::vector<SparseMatrix>& MM)
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
		fid >> m_InputData.Dimensions 
			>> m_InputData.NumberOfMassMtx
			>> m_InputData.NumberOfEigenvalues
			>> m_InputData.Tolerance
			>> m_InputData.MaxIter
			>> m_InputData.Sigma
			>> m_InputData.NumberOfKrylovBasis;

		// Set the matrices
		K0 = SparseMatrix(m_InputData.Dimensions, m_InputData.Dimensions);
		MM.reserve(m_InputData.NumberOfMassMtx);

		// Read the stiffness matrix K0
		for (int ii = 0; ii < m_InputData.Dimensions; ii++)
		{
			for (int jj = 0; jj < m_InputData.Dimensions; jj++)
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
		for (int im = 0; im < m_InputData.NumberOfMassMtx; im++)
		{
			SparseMatrix Mtemp(m_InputData.Dimensions, m_InputData.Dimensions);
			for (int ii = 0; ii < m_InputData.Dimensions; ii++)
			{
				for (int jj = 0; jj < m_InputData.Dimensions; jj++)
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
			MM.emplace_back(Mtemp);
		}

	}

	//Close file
	fid.close();
}

// Print info about the problem
void inverseFreeKrylovNLEigenSolver::printInfo()
{
	std::cout << "INFO" << std::endl;
	std::cout << "=====================================================" << std::endl;
	std::cout << "Dimension(s): " << m_InputData.Dimensions << std::endl;
	std::cout << "Number of mass matrices: " << m_InputData.NumberOfMassMtx << std::endl;
	std::cout << "Number of eigenpairs to compute: " << m_InputData.NumberOfEigenvalues << std::endl;
	std::cout << "Tolerance: " << m_InputData.Tolerance << std::endl;
	std::cout << "Max. iter: " << m_InputData.MaxIter << std::endl;
	std::cout << "Deflation value, sig: " << m_InputData.Sigma << std::endl;
	std::cout << "Number of Krylov basis, m: " << m_InputData.NumberOfKrylovBasis << std::endl;
	std::cout << "=====================================================\n\n";
	//ndim nm ne TOL max_iter sig nbasis

}

void inverseFreeKrylovNLEigenSolver::printResults(Vector& Omega, DenseMatrix& Phi) const
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

template<typename MtxType>
void inverseFreeKrylovNLEigenSolver::getFreqDependentStiffMtx(const MtxType& K0, const std::vector<MtxType>& MM, MtxType& Kn, data_type omega)
{
	//Initialize
	Kn = K0;

	for (int jj = 1; jj < m_InputData.NumberOfMassMtx; jj++)
	{
		Kn += (jj)*pow(omega, jj + 1.0) * MM[jj];
	}
}

template<typename MtxType>
void inverseFreeKrylovNLEigenSolver::getFreqDependentMassMtx(const std::vector<MtxType>& MM, MtxType& Mn, data_type omega)
{
	//Initialize
	Mn = MM[0];

	for (int jj = 1; jj < m_InputData.NumberOfMassMtx; jj++)
	{

		Mn += (jj + 1.0) * pow(omega, jj) * MM[jj];
	}
}

template<typename MtxType>
void inverseFreeKrylovNLEigenSolver::getGeneralizedFreqDependentMassMtx(const std::vector<MtxType>& MM, MtxType& Mlrls, data_type lr, data_type ls)
{
	//Initialize
	Mlrls.setZero();

	for (int jj = 0; jj < m_InputData.NumberOfMassMtx; jj++)
	{
		for (int kk = 0; kk < jj + 1; kk++)
		{
			Mlrls += pow(lr, kk) * pow(ls, jj - kk) * MM[jj];
		}
	}
}


void inverseFreeKrylovNLEigenSolver::generalizedDeflationProcedure(DenseMatrix& K0, std::vector<DenseMatrix>& MM_Update, DenseMatrix& Mlrls, data_type omega, Vector& phi)
{
	// Compute the generalized mass matrix Mlrls
	getGeneralizedFreqDependentMassMtx(MM_Update, Mlrls, omega, m_InputData.Sigma);

	// Check the number of mass matrices
	data_type aux = 0., sig = m_InputData.Sigma;
	// 1 MM
	if (m_InputData.NumberOfMassMtx == 1)
	{
		// For the case of 1 mm, only K0 is updated!
		aux = (sig - omega) / (phi.transpose() * Mlrls * phi);
		K0 += aux * MM_Update[0] * phi * phi.transpose() * MM_Update[0];

		return;
	}

	// 2 MM
	if (m_InputData.NumberOfMassMtx == 2)
	{
		// For the case of 2MM, all the terms 
		aux = (sig - omega) / (phi.transpose() * Mlrls * phi);
		auto phiphiT = phi * phi.transpose();

		// Error allocate dense to sparse in Eigen (To check and correct!!!!)
		// Update K0
		K0 += aux * (MM_Update[0] * phiphiT * MM_Update[0] + omega * (MM_Update[0] * phiphiT * MM_Update[1]
			+ MM_Update[1] * phiphiT * MM_Update[0]) + omega * omega * MM_Update[1] * phiphiT * MM_Update[1]);
		// Update M1
		MM_Update[0] -= aux * (MM_Update[0] * phiphiT * MM_Update[1]
			+ MM_Update[1] * phiphiT * MM_Update[0] + 2.0 * omega * MM_Update[1] * phiphiT * MM_Update[1]);
		// Update M2 
		MM_Update[1] -= aux * MM_Update[1] * phiphiT * MM_Update[1];

		return;
	}

}

// Compute the smallest eigenpair using modified free inverse Krylov method
void inverseFreeKrylovNLEigenSolver::computeSmallestEigenpairKrylov(const DenseMatrix& K0, const std::vector<DenseMatrix>& MM_Update, const std::vector<SparseMatrix>& MM_0, 
	data_type& omega, Vector& phi, int numberOfKrylovBasis, resultStatus& status)
{
	// Initialize
	DenseMatrix Ck(m_InputData.Dimensions, m_InputData.Dimensions);
	SparseMatrix Mlambda(m_InputData.Dimensions, m_InputData.Dimensions);
	DenseMatrix Am(numberOfKrylovBasis, numberOfKrylovBasis);
	DenseMatrix Zm(m_InputData.Dimensions, numberOfKrylovBasis);
	DenseMatrix Qm(numberOfKrylovBasis, numberOfKrylovBasis);
	Vector eVector(numberOfKrylovBasis); eVector.setOnes();
	double* diag = new double[numberOfKrylovBasis];
	double* subdiag = new double[numberOfKrylovBasis];
	data_type conv;
	double mu;

	getFreqDependentStiffMtx<DenseMatrix>(K0, MM_Update, Ck, omega);
	getFreqDependentMassMtx<SparseMatrix>(MM_0, Mlambda, omega);

	omega = RayleighQuotient<DenseMatrix>(Ck, Mlambda, phi);

	getFreqDependentStiffMtx<DenseMatrix>(K0, MM_Update, Ck, omega);
	getFreqDependentMassMtx<SparseMatrix>(MM_0, Mlambda, omega);

	Ck -= omega * Mlambda;

	LOG_INFO("In computeSmallestEigenpairKrylov");
	for (int k = 0; k < m_InputData.MaxIter; k++)
	{
		conv = phi.transpose() * phi;
		conv = sqrt(conv);
		// Construct a basis Zm using orthogonal basis by Arnoldi algorithm
		orthoArnoldiAlgorithm(Ck, Mlambda, numberOfKrylovBasis, phi, Zm);
		// Compute Am
		Am = Zm.transpose() * Ck * Zm;
		orthoLanczosAlgorithm(diag, subdiag, Am, numberOfKrylovBasis, eVector, Qm);

		// Find the smallest eigenpair of tridiagonal matrix
		eVector.setZero();
		smallestEigenpairTriDiagMtx(diag, subdiag, numberOfKrylovBasis, mu, eVector);

		// Update eigenpair solution
		omega += mu;
		phi = Zm * Qm * eVector;

		getFreqDependentStiffMtx<DenseMatrix>(K0, MM_Update, Ck, omega);
		getFreqDependentMassMtx<SparseMatrix>(MM_0, Mlambda, omega);
		Ck -= omega * Mlambda;
		
		// Check error criteria
		//conv = (Ck * phi).norm();
		
		data_type phiNorm = phi.transpose() * phi;
		conv = abs(1.0 - conv / (sqrt(phiNorm)));
		//LOG_INFO("Iter: {0}, conv = {1}", k, conv);
		if (conv < m_InputData.Tolerance)
		{
			phi *= 1.0/(phi.transpose() * Mlambda * phi);

			// Set status
			status.Convergence = conv;
			status.numberOfIterations = k;
			status.status = true;

			// Free memory
			delete[] diag;
			delete[] subdiag;
			return;
		}
			
	}

	// Free memory
	delete[] diag;
	delete[] subdiag;

	// Set status
	status.Convergence = (double)conv;
	status.numberOfIterations = m_InputData.MaxIter;
	status.status = false;

	LOG_WARN("conv = {0}", status.Convergence);
	LOG_ASSERT(false, "Error: It has reached the max. number of iterations!!");

}

// This functions compute the smallest eigenpairs of tridiagonal matrix 
// by the MRRR algorithm
void inverseFreeKrylovNLEigenSolver::smallestEigenpairTriDiagMtx(double* diag, double* subdiag, int numberOfBasis, double& eigValue, Vector& eigVector)
{
    // Checking the use of dstegr_
     /*dstegr_(JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
         *ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
         *LIWORK, INFO);*/
    integer il = 0, iu = 0, nin = 6, nout = 6, info;
    char jobz = 'V';
    char range = 'A';
    double abstol = 1e-14;
    double vl = 0.0, vu = 0.0;
    integer m = numberOfBasis, n = numberOfBasis;
    integer ldz, liwork, lwork;
    ldz = n;
    liwork = 10 * n;
    lwork = 18 * n;
    double* eigenvalues = new double[n];
    double* work = new double[lwork];
    double* eigenvector = new double[ldz * n];
    integer* isuppz = new integer[(2 * n)];
    integer* iwork = new integer[liwork];

    // Call DSTEGR subroutine (Algorithm to eigenlinear problem of symmetric tridiagonal  matrix)
    // MRRR algorithm
    dstegr_(&jobz, &range, &n, diag, subdiag, &vl, &vu, &il, &iu, &abstol, &m, eigenvalues, eigenvector, &ldz, isuppz, work,
        &lwork, iwork, &liwork, &info);

    // Get the smallest eigenpair
    eigValue = eigenvalues[0];
    for (int i = 0; i < n; i++)
    {
        eigVector(i) = eigenvector[i];
    }

    // free memory
    delete[] eigenvalues;
    delete[] work;
    delete[] eigenvector;
    delete[] isuppz;
    delete[] iwork;

}

template<typename MtxType>
data_type inverseFreeKrylovNLEigenSolver::RayleighQuotient(const MtxType& A, const MtxType& B, const Vector& eigVector)
{
    data_type result = 0.0;
    data_type num = eigVector.transpose() * A * eigVector;
    data_type den = eigVector.transpose() * B * eigVector;
    result = num / den;
	return result;
}

// Construct the Krylov basis by Lanczos procedure
void inverseFreeKrylovNLEigenSolver::orthoLanczosAlgorithm(double* diag, double* subdiag, const DenseMatrix& Ck, int numberOfBasis, const Vector& eigVector, DenseMatrix& Qm)
{
	// Get info and initialize	
	data_type norm = 0, alpha = 0, beta = 0.0;
	Vector  w(numberOfBasis); w.setZero();
	Qm.setZero();
	norm = eigVector.transpose() * eigVector;
	norm = sqrt(norm);
	Qm.col(0) = 1.0 / norm * eigVector;

	// Lanczos procedure
	for (int i = 0; i < numberOfBasis - 1; i++)
	{
		w = Ck * Qm.col(i);
		if (i > 0)
		{
			w -= beta * Qm.col(i - 1);
		}
		alpha = Qm.col(i).transpose() * w;
		w -= alpha * Qm.col(i);

		// Reorthogonalizing Qm (as suggested by Golub & Ye)
		if (numberOfBasis > 5)
		{
			for (int k = 0; k < i + 1; k++)
			{
				data_type temp = Qm.col(k).transpose() * w;
				w -= temp * Qm.col(k);
			}
		}

		beta = sqrt(w.transpose() * w);
		Qm.col(i + 1) = 1 / beta * w;

		// Set the values of the tridiagonal matrix
		diag[i] = (double)alpha;
		subdiag[i] = (double)beta;
	}

	// Computing the last component of the tridiagonal matrix
	w = Ck * Qm.col(numberOfBasis - 1);
	w -= beta * Qm.col(numberOfBasis - 2);
	data_type tempvalue = (Qm.col(numberOfBasis - 1).transpose() * w);
	diag[numberOfBasis - 1] = (double)tempvalue;
}

// Construct the Krylov basis by Arnoldi procedure
void inverseFreeKrylovNLEigenSolver::orthoArnoldiAlgorithm(const DenseMatrix& Ck, const DenseMatrix& B, int numberOfBasis, const Vector& eigVector, DenseMatrix& Zm)
{
	// Get info and initialize
	Vector w(numberOfBasis);
	DenseMatrix Hm(numberOfBasis, numberOfBasis);
	data_type norm = 0;

	norm = sqrt(eigVector.transpose() * B * eigVector);
	Zm.col(0) = 1 / norm * eigVector;

	// B-orthonormal basis by the Arnoldi algorithm
	for (int i = 0; i < numberOfBasis - 1; i++)
	{
		w = Ck * Zm.col(i);
		for (int j = 0; j < i + 1; j++)
		{
			Hm(j, i) = Zm.col(j).transpose() * B * w;
			w -= Hm(j, i) * Zm.col(j);
		}
		norm = sqrt(w.transpose() * B * w);
		Zm.col(i + 1) = 1 / norm * w;
	}

}
