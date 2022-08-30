#pragma once
#include "NLEigenSolver.h"

//  Implementation with support to sparse matrix (Optmizing the code!!!!)
class JacobiDavidsonNLEigenSolver : public NLEigenSolver
{
public:
	JacobiDavidsonNLEigenSolver(const std::string& filepath);
	~JacobiDavidsonNLEigenSolver();

	virtual bool execute() override;
	//Miscellaneous
	virtual bool findEigenvaluesFromInitialGuess() override;

private:
	void readFileAndGetStiffMassMatrices(SparseMatrix& K0, std::vector<SparseMatrix>& MM, Vector& Omega);
	void printResults(Vector& Omega, DenseMatrix& Phi) const;
	void getFreqDependentStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, SparseMatrix& Kn, data_type omega);  // Kn(lr)
	void getFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mn, data_type omega);                               // Mn(lr)
	void getGeneralizedFreqDependentMassMtx(const std::vector<SparseMatrix>& MM, SparseMatrix& Mlrls, data_type lr, data_type ls);         // M(lr,ls)
	void getEffectiveStiffMtx(const SparseMatrix& K0, const std::vector<SparseMatrix>& MM, DenseMatrix& Keff, data_type omega);     // Keff = K(lr)-lr*M(lr)
	void projectEffectiveStiffMatrix(DenseMatrix& Keff, DenseMatrix& B_s, int indexEig);
	bool iterativeLinearSolver(const DenseMatrix& A, const Vector& b, Vector& x);
	void UpdateVectorSolution(DenseMatrix& Keff, DenseMatrix& Phi, DenseMatrix& B_r, int index);


private:
	int m_Dimensions;
	int m_NumberOfMassMtx;
	int m_NumberOfEigenValues;
	int m_MaxIter;
	double m_TOL;
	std::string m_FilePath;
	bool m_hasInitialTrial = false;

};


