#pragma once 

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <iomanip>

#include "Log.h"

//#define ENABLE_BLAZE

#ifdef ENABLE_BLAZE
	#include <blaze/Math.h>
#else
	#include <Eigen/Dense>
	#include <Eigen/Core>
#endif




