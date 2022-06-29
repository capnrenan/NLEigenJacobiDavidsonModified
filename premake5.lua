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
		"%{prj.name}/vendor/mpfr/include/**.h",
		"%{prj.name}/vendor/mpir/include/**.h"
	}
	
	
	pchheader "nlpch.h"
	pchsource "NLEigenSolver/src/nlpch.cpp"
	

	-- Additional library directories
	libdirs { "%{prj.name}/vendor/CLAPACK/lib",
				"%{prj.name}/vendor/mpfr/lib",
				"%{prj.name}/vendor/mpir/lib"}
	

				
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
		"%{prj.name}/vendor/mpfr/include",
		"%{prj.name}/vendor/mpir/include"
	}

	-- Filter: Configurations only applied to specific platforms
	filter "system:windows"
		systemversion "latest"
		
	filter "configurations:Debug"
		defines { "_DEBUG" }
		runtime "Debug"
		symbols "On"
		links { "libf2cd.lib","blasd.lib","lapackd.lib", "mpir.lib", "mpfr.lib"}
	  
	filter "configurations:Release"
		defines { "_RELEASE" }
		runtime "Release"
		optimize "On"
		symbols "On"
		
		links { "libf2c.lib","blas.lib","lapack.lib","mpir.lib", "mpfr.lib"}	  

