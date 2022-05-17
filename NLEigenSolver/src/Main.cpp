#include "nlpch.h"
#include "NLEigenSolver.h"

// Entry Point
int main(int argc, char* argv[])
{
	PROFILE_BEGIN_SESSION("Eigenvalue");
	std::shared_ptr<NLEigenSolver> app = NLEigenSolver::Create(argv[1]);
	bool status = app->execute();
	PROFILE_END_SESSION();

	// Check status
	if (!status)
		return 1;

	return 0;
}