#include "nlpch.h"
#include "NLEigenJacobiDavidson.h"


int main(int argc, char* argv[])
{
	NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson("examples/K.dat");
	app->execute();

	delete app;
	return 0;
}