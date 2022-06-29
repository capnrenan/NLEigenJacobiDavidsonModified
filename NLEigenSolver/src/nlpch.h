#pragma once 

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <iomanip>

#include "Log.h"
#include "Instrumentor.h"

#include <mpreal.h>



//#define ENABLE_BLAZE

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




