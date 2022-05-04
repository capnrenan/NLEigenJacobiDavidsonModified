#pragma once
#include "NLEigenSolver.h"

#include <blaze/Math.h>

class BlazeNLEigenSolver : public NLEigenSolver
{
public:
	BlazeNLEigenSolver(const std::string& filepath);
	virtual ~BlazeNLEigenSolver();

	virtual bool execute() override;

private:
	void readFileAndGetStiffMassMatrices(blaze::DynamicMatrix<double, blaze::rowMajor>& K0, std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM);
	void printResults(blaze::DynamicVector<double, blaze::columnVector>& Omega, blaze::DynamicMatrix<double, blaze::rowMajor>& Phi) const;

	void getFreqDependentStiffMtx(const blaze::DynamicMatrix<double, blaze::rowMajor>& K0, const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Kn, double omega);  // Kn(lr)
	void getFreqDependentMassMtx(const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Mn, double omega);                               // Mn(lr)
	void getGeneralizedFreqDependentMassMtx(const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Mlrls, double lr, double ls);         // M(lr,ls)
	void getEffectiveStiffMtx(const blaze::DynamicMatrix<double, blaze::rowMajor>& K0, const std::vector<blaze::DynamicMatrix<double, blaze::rowMajor>>& MM, blaze::DynamicMatrix<double, blaze::rowMajor>& Keff, double omega);     // Keff = K(lr)-lr*M(lr)
	void projectEffectiveStiffMatrix(blaze::DynamicMatrix<double, blaze::rowMajor>& Keff, blaze::DynamicMatrix<double, blaze::rowMajor>& B_s, int indexEig);
	bool iterativeLinearSolver(const blaze::DynamicMatrix<double, blaze::rowMajor>& A, const blaze::DynamicVector<double, blaze::columnVector>& b, blaze::DynamicVector<double, blaze::columnVector>& x);

public:
	int m_Dimensions;
	int m_NumberOfMassMtx;
	int m_NumberOfEigenValues;
	int m_MaxIter;
	double m_TOL;
	std::string m_FilePath;
};