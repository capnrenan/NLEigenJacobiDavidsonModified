#pragma once
#include "nlpch.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>


#define QUAD_PRECISION 0
// Check the quad precision
#if QUAD_PRECISION
#include <unsupported/Eigen/MPRealSupport>

// Set precision (long double - 128 bits)
const int digits = 20;
mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
using data_type = mpfr::mpreal;
using DenseMatrix = Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Vector<mpfr::mpreal, Eigen::Dynamic>;

using SparseMatrix = Eigen::SparseMatrix<mpfr::mpreal>;

#else
using data_type = double;
using DenseMatrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

using SparseMatrix = Eigen::SparseMatrix<data_type>;

#endif


#define PROFILING 0
#if PROFILING
#define PROFILE_BEGIN_SESSION(...) Instrumentor::Get().BeginSession(__VA_ARGS__)
#define PROFILE_END_SESSION() Instrumentor::Get().EndSession()
#define PROFILE_SCOPE(name) IntrumentationTimer timer##__LINE__(name)
//#define PROFILE_FUNCTION() PROFILE_SCOPE(__FUNCTION__)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__FUNCSIG__)
#else
#define PROFILE_BEGIN_SESSION(...) 
#define PROFILE_END_SESSION() 
#define PROFILE_SCOPE(name)
//#define PROFILE_FUNCTION() PROFILE_SCOPE(__FUNCTION__)
#define PROFILE_FUNCTION() 
#endif


class NLEigenMethods
{
public:
	enum class Method
	{
		None = 0, JacobiDavidson = 1, inverseFreeKrylov = 2
	};

};

class NLEigenSolver
{
public:
	virtual ~NLEigenSolver() {};

	virtual bool execute() = 0;
	virtual bool findEigenvaluesFromInitialGuess() = 0;

	static std::shared_ptr<NLEigenSolver>Create(const NLEigenMethods::Method& method, const std::string& filepath);
};