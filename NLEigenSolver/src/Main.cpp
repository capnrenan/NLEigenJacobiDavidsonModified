#include "nlpch.h"
#include "NLEigenSolver.h"


// Entry Point
int main(int argc, char* argv[])
{
	// Set precision (long double - 128 bits)
	mpfr::mpreal::set_default_prec(128);

	PROFILE_BEGIN_SESSION("Eigenvalue routine");
	std::shared_ptr<NLEigenSolver> app = NLEigenSolver::Create(argv[1]);
	bool status = app->execute();
	PROFILE_END_SESSION();

	// Check status
	if (!status)
		return 1;

	return 0;
}