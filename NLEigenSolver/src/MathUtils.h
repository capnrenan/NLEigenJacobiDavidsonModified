#include "nlpch.h"
#include<functional>


namespace MathUtils
{
	// template Linear regula-falsi function
	// DataFormat - double, float, mpfr::mpreal;
	
	template<typename DataFormat>
	DataFormat linearRegularFalsi(DataFormat initialStep, const std::function<DataFormat(DataFormat)>& function, 
		bool& status, double TOL = 1e-8, uint32_t MAX_ITER = 20)
	{
		// Initializing
		status = true;
		DataFormat x0, f0, x1, f1, delta, result, fCurrent,	 residual;
		x0 = initialStep; 
		f0 = function(x0); 
		delta = abs(x0) / 5.0;
		x1 = x0 - delta; f1 = function(x1);
		residual = 1.0;
		result = x0;
		
		// Computing the root
		for (uint32_t i = 0; i < MAX_ITER; i++)
		{
			DataFormat dF = f0 - f1;
			// Check dF 
			if (abs(dF) < 1e-12)
			{
				LOG_ERROR( "No local solution, dF = {0}", dF);
				status = false;
				return result;
			}
			else
			{
				result = (f0 * x1 - f1 * x0) / dF;
			}

			// Compute the residual
			// Note: The search algorithm does not work with multiple roots
			fCurrent = function(result);
			residual = abs(result - x0) / result;
			residual = (abs((result - x1) / result) < residual) ? abs((result - x1) / result) : residual;
			//residual = (abs(fCurrent) < residual) ? abs(fCurrent) : residual;
			//residual = abs(fCurrent);
			LOG_INFO("iter: {0}, error = {1}, value = {2}", i, residual, result);

			// Check residual
			if (residual < TOL)
				return result;

			// Update x0, x1, f0, f1
			if (result < x0)
			{
				x1 = x0;
				f1 = f0;
				x0 = result;
				f0 = fCurrent;
			}
			else if (result < x1)
			{
				x1 = result;
				f1 = fCurrent;
			}
			else
			{
				x0 = x1;
				f0 = f1;
				x1 = result;
				f1 = fCurrent;
			}
		}
		
		LOG_WARN("Algorithm has reach the maximum number of iterations!!!! ");
		status = false;
		return result;
	}
}