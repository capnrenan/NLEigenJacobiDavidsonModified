#include "nlpch.h"
#include "NLEigenJacobiDavidson.h"


int main(int argc, char* argv[])
{
	//std::string filepath;
	//std::cin >> filepath;
	//NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson("examples/K.dat");
	NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson(argv[1]);
	bool status = app->execute();

	delete app;
	if (status)
		return 0;
	else
		return 1;
}