#include "nlpch.h"
#include "NLEigenSolver.h"


// Entry Point
int main(int argc, char* argv[])
{
	//NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson(argv[1]);
	//std::shared_ptr<NLEigenSolver> app = NLEigenSolver::Create("examples/K.txt");
	Instrumentor::Get().BeginSession("Profile");
	//auto app = NLEigenSolver::Create(argv[1]);
	std::shared_ptr<NLEigenSolver> app = NLEigenSolver::Create("examples/K.txt");
	bool status = app->execute();
	Instrumentor::Get().EndSession();

	if (status)
		return 0;
	else
		return 1;

}