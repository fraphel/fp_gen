#ifndef __CardioXcomp__output__
#define __CardioXcomp__output__

#include "platform.h"
#include <iostream>
#include <iomanip>
#include <sstream>

#include <vector>
#include <petsc.h>
#include <fstream>

#include "mesh.hpp"
#include "dof.hpp"

using namespace std;

namespace cardioxcomp{
  
  void writeSolution_0d_Headers(std::vector<string> name, string filename,bool binPostProc=false);
  void writeSolution_0d(std::vector<Vec> vecvalue, string filename, double time,bool binPostProc=false);
  void writeSolution_0d(double time, double value, string filename,bool binPostProc=false);


  void InitEnsight(const std::string& filename, Mesh* mesh, std::vector<string> varnames, double tMax, double dt, int nMax, int freqWrite);
  void InitEnsight(const std::string& filegeoname, const std::string& filename, Mesh* mesh, std::vector<string> varnames, double tMax, double dt, int nMax, int freqWrite);

  void writeEnsight(string filename, int numIt,double time, Vec vecSol, std::vector<string> varnames, Dof* dof);
  void writeEnsight(string filename, int numIt,double time, std::vector<Vec> vecSol, std::vector<string> varnames, Dof* dof);
  
}


#endif
