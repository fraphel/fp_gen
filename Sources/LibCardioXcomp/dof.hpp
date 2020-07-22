#ifndef __CardioXcomp__dof__
#define __CardioXcomp__dof__

#include <vector>
#include <petscsys.h>
#include <petscao.h>
#include <petscmat.h>
#include <parmetis.h>

#include "platform.h"
#include "mesh.hpp"
#include "csrmatrixpattern.hpp"

namespace cardioxcomp {

  class Dof {
    const Mesh& m_mesh;
    //! number of components per vertex
    int m_numComp;
    //! number of dof = number of vertices times m_numComp
    int m_numDof;
    //! number of dof on this proc
    int m_numDofLocal;
    //! CSR matrix pattern
    CSRMatrixPattern m_pattern;
    //! m_dofPart[idof] = the rank of the processor which keeps in memory the dof idof
    std::vector<int> m_dofPart;
    //! m_eltPart[iele] = rank of the processor which keep in memory the element iele
    std::vector<int> m_eltPart;
    //! mapping between global problem numbering and matrix dof numbering
    AO m_ao;
    bool m_initAO;
    //! mapping between local to global ordering of dof.
    ISLocalToGlobalMapping m_mappingNodes;
    bool m_initMappingNodes;
  public:
    Dof(const Mesh& mesh,int numComp=1);
    ~Dof();
    void initializePattern(int sizeProc, int rankProc);
    void buildPattern(int rankProc, const PetscInt* dofRepartition);
    void partitionDof(idx_t numPart,int rankPart,MPI_Comm comm);
    const CSRMatrixPattern& pattern(){return m_pattern;}
    const std::vector<int>& eltPart(){return m_eltPart;}
    void allocateMatrixVec(int size, int rank,Mat* matrix,int numMat,Vec* vec,int numVec);
    AO& ao(){return m_ao;}
    const AO& ao() const {return m_ao;}
    std::vector<int>& dofPart(){return m_dofPart;}
    const std::vector<int>& dofPart() const {return m_dofPart;}
    const int& numComp() const {return m_numComp;}
    const int& numDof() const {return m_numDof;}
    const int& numDofLocal() const {return m_numDofLocal;}

#ifdef PETSC_3_6
    const ISLocalToGlobalMapping& locToGlob() const {return m_mappingNodes;}
#else
    ISLocalToGlobalMapping return_mappingNodes() const {return m_mappingNodes;}
#endif

  };
}

#endif

