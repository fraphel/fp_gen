#pragma once

#include <petscsys.h>
#include "bidomainproblem.hpp"
#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>

#ifdef _WIN32
#include "problem.hpp"
#include "luafile.hpp"
#else
#include "../../LibCardioXcomp/problem.hpp"
#include "../../LibCardioXcomp/luafile.hpp"
#endif
//#include "run0D.h"

using namespace std;

namespace cardioxcomp {

#ifdef _WIN32
	extern "C" __declspec(dllexport)
#endif
 void initializeCardioXComp(char*,int&,int&);//configfile,nb el, nbit

#ifdef _WIN32
	extern "C" __declspec(dllexport)
#endif
  void runCardioXCompFiles(char*,char*,double*);//configfile, variablesfile, solution

#ifdef _WIN32
	extern "C" __declspec(dllexport)
#endif
  void runCardioXCompUI(LuaFile*,LuaVariables*,double*);

#ifdef _WIN32
	extern "C" __declspec(dllexport)
#endif
  Problem* runCardioXCompUI_First(LuaFile*,LuaVariables*,double*, double, int);//From 0 to tmax

#ifdef _WIN32
	extern "C" __declspec(dllexport)
#endif
  Problem* runCardioXCompUI_Others(double*, int, double, Problem*, int, int& sizeData);// from tlast to tnext

#ifdef _WIN32
	extern "C" __declspec(dllexport)
#endif
  Problem* runCardioXCompUI_FirstCpp(LuaFile*,LuaVariables*,double*, double, int, int, Problem*);//From 0 to tmax

}
