#include "nlpch.h"
#include "NLEigenJacobiDavidson.h"

// Entry Point
int main(int argc, char* argv[])
{
	NLEigenJacobiDavidson* app = new NLEigenJacobiDavidson(argv[1]);
	bool status = app->execute();

	delete app;
	if (status)
		return 0;
	else
		return 1;

}