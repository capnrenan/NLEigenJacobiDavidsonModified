#include "nlpch.h"
#include "Log.h"

std::shared_ptr<spdlog::logger> Log::s_CoreLogger;

void Log::Init()
{
	//[%T] - timestamp
	//%n - name of the Logger
	//%v% - message
	spdlog::set_pattern("%^[%T] %n: %v%$");
	s_CoreLogger = spdlog::stdout_color_mt("NLEIGEN");
	s_CoreLogger->set_level(spdlog::level::trace);

}

