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
	LOG_INFO("Reading filedata...\n");

	SparseMatrix K0;
	std::vector<SparseMatrix> MM;
	Vector Omega;
	readFileAndGetStiffMassMatrices(K0, MM, Omega);

	

	return false;
}

bool inverseFreeKrylovNLEigenSolver::findEigenvaluesFromInitialGuess()
{
	return false;
}

void inverseFreeKrylovNLEigenSolver::readFileAndGetStiffMassMatrices(SparseMatrix& K0, std::vector<SparseMatrix>& MM, Vector& Omega)
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
		fid >> m_InputData->Dimensions >> m_InputData->NumberOfMassMtx  >> m_InputData->NumberOfEigenvalues >> m_InputData->Tolerance;

		// Set the matrices
		K0 = SparseMatrix(m_InputData->Dimensions, m_InputData->Dimensions);
		Omega = Vector(m_InputData->NumberOfEigenvalues);
		Omega.setZero();
		MM.reserve(m_InputData->NumberOfMassMtx);

		// Read the stiffness matrix K0
		for (int ii = 0; ii < m_InputData->Dimensions; ii++)
		{
			for (int jj = 0; jj < m_InputData->Dimensions; jj++)
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
		for (int im = 0; im < m_InputData->NumberOfMassMtx; im++)
		{
			SparseMatrix Mtemp(m_InputData->Dimensions, m_InputData->Dimensions);
			for (int ii = 0; ii < m_InputData->Dimensions; ii++)
			{
				for (int jj = 0; jj < m_InputData->Dimensions; jj++)
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

	}

	//Close file
	fid.close();
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

void inverseFreeKrylovNLEigenSolver::getFreqDependentStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, SparseMatrix& Kn, data_type omega)
{
	//Initialize
	Kn = K0;

	for (int jj = 1; jj < m_InputData->NumberOfMassMtx; jj++)
	{
		Kn += (jj)*pow(omega, jj + 1.0) * MM[jj];
	}
}

void inverseFreeKrylovNLEigenSolver::getFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mn, data_type omega)
{
	//Initialize
	Mn = MM[0];

	for (int jj = 1; jj < m_InputData->NumberOfMassMtx; jj++)
	{

		Mn += (jj + 1.0) * pow(omega, jj) * MM[jj];
	}
}

void inverseFreeKrylovNLEigenSolver::getGeneralizedFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mlrls, data_type lr, data_type ls)
{
	//Initialize
	Mlrls.setZero();

	for (int jj = 0; jj < m_InputData->NumberOfMassMtx; jj++)
	{
		for (int kk = 0; kk < jj + 1; kk++)
		{
			Mlrls += pow(lr, kk) * pow(ls, jj - kk) * MM[jj];
		}
	}
}

// Compute the smallest eigenpair using modified free inverse Krylov method
void inverseFreeKrylovNLEigenSolver::computeSmallestEigenpairKrylov(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM_Update, const std::vector<SparseMatrix>& MM_0, 
	double& omega, Vector& phi, int numberOfKrylovBasis)
{
	// Initialize
	SparseMatrix Ck(m_InputData->Dimensions, m_InputData->Dimensions);
	SparseMatrix Mlambda(m_InputData->Dimensions, m_InputData->Dimensions);
	SparseMatrix Am(numberOfKrylovBasis, numberOfKrylovBasis);
	DenseMatrix Zm(m_InputData->Dimensions, numberOfKrylovBasis);
	DenseMatrix Qm(numberOfKrylovBasis, numberOfKrylovBasis);
	Vector eVector(numberOfKrylovBasis); eVector.setOnes();
	double* diag = new double[numberOfKrylovBasis];
	double* subdiag = new double[numberOfKrylovBasis];
	double mu, conv;

	getFreqDependentStiffMtx(K0, MM_Update, Ck, omega);
	getFreqDependentMassMtx(MM_0, Mlambda, omega);
	Ck -= omega * Mlambda;
	omega = RayleighQuotient(Ck, Mlambda, phi);

	for (int k = 0; k < m_InputData->MaxIter; k++)
	{
		// Construct a basis Zm using orthogonal basis by Arnoldi algorithm
		orthoArnoldiAlgorithm(Ck, Mlambda, numberOfKrylovBasis, phi, Zm);
		// Compute Am
		Am = Zm.transpose() * Ck * Zm;
		orthoLanczosAlgorithm(diag, subdiag, Ck, numberOfKrylovBasis, eVector, Qm);

		// Find the smallest eigenpair of tridiagonal matrix
		smallestEigenpairTriDiagMtx(diag, subdiag, numberOfKrylovBasis, mu, eVector);

		// Update eigenpair solution
		omega += mu;
		phi = Zm * Qm * eVector;

		getFreqDependentStiffMtx(K0, MM_Update, Ck, omega);
		getFreqDependentMassMtx(MM_0, Mlambda, omega);
		Ck -= omega * Mlambda;
		
		// Check error criteria
		conv = (Ck * phi).norm();
		if (conv < m_InputData->Tolerance)
		{
			phi *= phi.transpose() * Mlambda * phi;
			return;
		}
			
	}

}

// This functions compute the smallest eigenpairs of tridiagonal matrix 
// by the MRRR algorithm
void inverseFreeKrylovNLEigenSolver::smallestEigenpairTriDiagMtx(double* diag, double* subdiag, int numberOfBasis, double& eigValue, Vector& eigVector)
{
    // Checking the use of dstegr_
    // To do: 
     /*dstegr_(JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
         *ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
         *LIWORK, INFO);*/
    integer il = 0, iu = 0, nin = 6, nout = 6, info;
    char jobz = 'V';
    char range = 'I';
    double abstol = 1e-10;
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
    for (uint32_t i = 0; i < n; i++)
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

data_type inverseFreeKrylovNLEigenSolver::RayleighQuotient(const SparseMatrix& A, const SparseMatrix& B, const Vector& eigVector)
{
    data_type result = 0.0;
    data_type num = eigVector.transpose() * A * eigVector;
    data_type den = eigVector.transpose() * B * eigVector;
    result = num / den;
	return result;
}

// Construct the Krylov basis by Lanczos procedure
void inverseFreeKrylovNLEigenSolver::orthoLanczosAlgorithm(double* diag, double* subdiag, const SparseMatrix& Ck, int numberOfBasis, const Vector& eigVector, DenseMatrix& Qm)
{
	// Get info and initialize	
	double norm = 0, alpha = 0, beta = 0.0;
	Vector temp(numberOfBasis), w(numberOfBasis);
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
		if (numberOfBasis > 6)
		{
			for (int k = 0; k < i + 1; k++)
			{
				temp = Qm.col(k).transpose() * w;
				w -= temp * Qm.col(k);
			}
		}

		beta = sqrt(w.transpose() * w);
		Qm.col(i + 1) = 1 / beta * w;

		// Set the values of the tridiagonal matrix
		diag[i] = alpha;
		subdiag[i] = beta;
	}

	// Computing the last component of the tridiagonal matrix
	w = Ck * Qm.col(numberOfBasis - 1);
	w -= beta * Qm.col(numberOfBasis - 2);
	diag[numberOfBasis - 1] = (Qm.col(numberOfBasis - 1).transpose() * w);
}

// Construct the Krylov basis by Arnoldi procedure
void inverseFreeKrylovNLEigenSolver::orthoArnoldiAlgorithm(const SparseMatrix& Ck, const SparseMatrix& B, int numberOfBasis, const Vector& eigVector, DenseMatrix& Zm)
{
	// Get info and initialize
	Vector w(numberOfBasis);
	DenseMatrix Hm(numberOfBasis, numberOfBasis);
	double norm = 0;

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
