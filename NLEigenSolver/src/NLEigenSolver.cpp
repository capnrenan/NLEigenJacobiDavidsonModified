#include "nlpch.h"
#include "NLEigenSolver.h"
#include "JacobiDavidsonNLEigenSolver.h"
#include "inverseFreeKrylovNLEigenSolver.h"

std::shared_ptr<NLEigenSolver> NLEigenSolver::Create(const NLEigenMethods::Method& method, const std::string& filepath)
{
	Log::Init();
	// Check the math library API
	switch (method)
	{
		case NLEigenMethods::Method::None:
		{
			LOG_ASSERT(false, "NLEigenMethods::Method::None is currently not supported!");
			return nullptr;
		}

		case NLEigenMethods::Method::JacobiDavidson:
		{
			return std::make_shared<JacobiDavidsonNLEigenSolver>(filepath);
		}

		case NLEigenMethods::Method::inverseFreeKrylov:
		{
			return std::make_shared<inverseFreeKrylovNLEigenSolver>(filepath);
		}
	}

	LOG_ASSERT(false, "Unknown NLEigensolver method!!!");
	return nullptr;

}
