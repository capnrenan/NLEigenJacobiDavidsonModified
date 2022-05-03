#include "nlpch.h"
#include "NLEigenSolver.h"
#include "BlazeNLEigenSolver.h"
#include "EigenNLEigenSolver.h"

class MathLibrary
{
public:
	enum class API
	{
		None = 0, Eigen = 1, Blaze = 2
	};

	inline static API GetAPI() { return s_API; };

private:
	static API s_API;
};

// Check if the math library API
#ifdef ENABLE_BLAZE
	MathLibrary::API MathLibrary::s_API = MathLibrary::API::Blaze;
#else
	MathLibrary::API MathLibrary::s_API = MathLibrary::API::Eigen;
#endif 


std::shared_ptr<NLEigenSolver> NLEigenSolver::Create(const std::string& filepath)
{
	Log::Init();
	// Check the math library API
	switch (MathLibrary::GetAPI())
	{
		case MathLibrary::API::None:
		{
			LOG_ASSERT(false, "MathLibrary::API::None is currently not supported!");
			return nullptr;
		}

		case MathLibrary::API::Eigen:
		{
			return std::make_shared<EigenNLEigenSolver>(filepath);
		}

		case MathLibrary::API::Blaze:
		{
			return std::make_shared<BlazeNLEigenSolver>(filepath);
		}
	}

	LOG_ASSERT(false, "Unknown Math Library!!!");
	return nullptr;

}
