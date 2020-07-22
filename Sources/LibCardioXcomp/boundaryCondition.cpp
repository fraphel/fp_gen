#include "stdafx.h"

#include <iostream>
#include "boundaryCondition.hpp"

#undef VERBOSE_DEBUG
#define TGV 1e20

namespace cardioxcomp
{
  
  void applyDirichletBC(Mat& mat,Vec& rhs,double val,int label,Mesh& mesh,const Dof& dof,int rankProc,int iComp){
#ifdef VERBOSE_DEBUG
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n========  Apply BC on label %d, component %d, on proc #%d ========\n",label,iComp,rankProc);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"List of vertices with label %d located on proc %d: ",label,rankProc);
    for (int i=0; i<mesh.vertexPerLabel[label].size(); i++) {
      int ivertex = mesh.vertexPerLabel[label][i];
      if(dof.dofPart()[ivertex] == rankProc){
        PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d ", mesh.vertexPerLabel[label][i]);
      }
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n.........  End Apply BC on label %d, component %d, on proc #%d .......\n",label,iComp,rankProc);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif

    for (size_t i=0; i<mesh.vertexPerLabel[label].size(); i++) {
      int ivertex = mesh.vertexPerLabel[label][i];
      platformInt idof = ivertex * dof.numComp() + iComp;

      AOApplicationToPetsc(dof.ao(), 1, &idof);
      MatSetValue(mat, idof, idof, TGV, ADD_VALUES);
      VecSetValue(rhs, idof, val*TGV, ADD_VALUES);

    }
  }
}
