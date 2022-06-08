projectName = "NLEigenSolver"
workspace "NLEigenSolver"
	architecture "x64"
	startproject  "NLEigenSolver"
	
	configurations
	{
		"Debug",
		"Release",
	}

outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"


-- Include directories relative to root folder (sol. directory)
IncludeDir = {}
--IncludeDir["Eigen"] = projectName .. "/vendor/Eigen"
--IncludeDir["Blaze"] = projectName .. "/vendor/Blaze"
--IncludeDir["CLAPACK"] = projectName .. "/vendor/CLAPACK/include"
--IncludeDir["spdlog"] = projectName .. "/vendor/spdlog/include"

project "NLEigenSolver"
	location "NLEigenSolver"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++17"
	staticruntime "off"

	targetdir ("bin/" .. outputdir .. "/%{prj.name}")
	objdir ("bin-int/" .. outputdir .. "/%{prj.name}")
	
	
	files
	{
		"%{prj.name}/src/**.h",
		"%{prj.name}/src/**.cpp",
	}
	
	
	pchheader "nlpch.h"
	pchsource "NLEigenSolver/src/nlpch.cpp"
	

	-- Additional library directories
	libdirs { "%{prj.name}/vendor/CLAPACK/lib" }
	
	-- OpenMP support
	buildoptions { "/openmp" }
	
	defines 
	{
		"NDEBUG", "_CONSOLE", 
	}
	
	includedirs
	{
		"%{prj.name}/src",
		"%{prj.name}/vendor/Eigen",
		"%{prj.name}/vendor/Blaze",
		"%{prj.name}/vendor/CLAPACK/include",
		"%{prj.name}/vendor/spdlog/include",
	}

	-- Filter: Configurations only applied to specific platforms
	filter "system:windows"
		cppdialect "C++17"
		systemversion "latest"
		
	filter "configurations:Debug"
      defines { "DEBUG" }
	  runtime "Debug"
      symbols "On"
	  links { "libf2cd.lib","blasd.lib","lapackd.lib"}
	  
	filter "configurations:Release"
      defines { "RELEASE" }
	  runtime "Release"
	  optimize "On"
      symbols "On"
	  links { "libf2c.lib","blas.lib","lapack.lib"}	  

