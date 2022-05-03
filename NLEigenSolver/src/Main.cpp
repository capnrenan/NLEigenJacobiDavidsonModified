#include "nlpch.h"
#include "NLEigenSolver.h"

// Entry Point
int main(int argc, char* argv[])
{
	//NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson(argv[1]);
	//std::shared_ptr<NLEigenSolver> app = NLEigenSolver::Create("examples/K.txt");
	auto app = NLEigenSolver::Create(argv[1]);
	bool status = app->execute();

	if (status)
		return 0;
	else
		return 1;

}