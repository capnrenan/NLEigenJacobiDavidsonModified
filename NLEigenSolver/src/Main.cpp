#include "nlpch.h"
#include "NLEigenJacobiDavidson.h"

// Entry Point
int main(int argc, char* argv[])
{
	NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson(argv[1]);
	//NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson("C:\MyPrograms\NLEigenSolver\NLEigenSolver\examples\K.txt");
	//NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson("examples/K.txt");
	bool status = app->execute();

	delete app;
	if (status)
		return 0;
	else
		return 1;

}