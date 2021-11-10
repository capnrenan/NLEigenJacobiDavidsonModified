#pragma once

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/fmt/ostr.h"


class Log
{
private:
	static std::shared_ptr<spdlog::logger> s_CoreLogger;
public:
	static void Init();

	inline static std::shared_ptr<spdlog::logger>& GetCoreLogger() { return s_CoreLogger; }

};

//Core log macros
#define LOG_ERROR(...)  ::Log::GetCoreLogger()->error(__VA_ARGS__)
#define LOG_WARN(...)   ::Log::GetCoreLogger()->warn(__VA_ARGS__)
#define LOG_INFO(...)   ::Log::GetCoreLogger()->info(__VA_ARGS__)
#define LOG_TRACE(...)  ::Log::GetCoreLogger()->trace(__VA_ARGS__)
#define LOG_FATAL(...)  ::Log::GetCoreLogger()->fatal(__VA_ARGS__)

#ifdef _DEBUG
	#define ENABLE_ASSERTS
#endif

#ifdef ENABLE_ASSERTS
#define LOG_ASSERT(x, ...) {if(!(x)) {LOG_ERROR(" Assertion Failed: {0}", __VA_ARGS__); __debugbreak();}}
#else
#define LOG_ASSERT(x, ...)
#endif