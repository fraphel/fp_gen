#ifndef __CardioXcomp__boundaryconditions__
#define __CardioXcomp__boundaryconditions__

#include <petscsys.h>
#include <petscmat.h>
#include <petscvec.h>

#include "mesh.hpp"
#include "dof.hpp"
#include "platform.h"

namespace cardioxcomp {

  void applyDirichletBC(Mat& mat,Vec& rhs,double val,int label,Mesh& mesh,const Dof& dof,int rankProc,int iComp=0);
  
}

#endif

