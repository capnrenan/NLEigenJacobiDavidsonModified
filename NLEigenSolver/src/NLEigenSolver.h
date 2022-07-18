#pragma once

#include <memory>
#include <string>

class NLEigenSolver
{
public:
	virtual ~NLEigenSolver() {};

	virtual bool execute() = 0;
	virtual bool findEigenvaluesFromInitialGuess() = 0;

	static std::shared_ptr<NLEigenSolver> Create(const std::string& filepath);
};