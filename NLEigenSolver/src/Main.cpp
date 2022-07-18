#include "nlpch.h"
#include "NLEigenSolver.h"



// Entry Point
int main(int argc, char* argv[])
{
	// Set precision (long double - 128 bits)
	const int digits = 30;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
	PROFILE_BEGIN_SESSION("Eigenvalue routine");
	std::shared_ptr<NLEigenSolver> app = NLEigenSolver::Create(argv[1]);
	bool status = app->execute();
	//bool status = app->findEigenvaluesFromInitialGuess();
	PROFILE_END_SESSION();

	// Check status
	if (!status)
		return 1;

	return 0;
}