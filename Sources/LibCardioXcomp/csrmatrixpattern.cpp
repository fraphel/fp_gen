#include "stdafx.h"

#include "csrmatrixpattern.hpp"

namespace cardioxcomp {
  std::ostream& operator <<  (std::ostream& f, const CSRMatrixPattern & pattern )  {
    f << "\n============ CSR Pattern ===========\n";
    f << " Number of rows = " << pattern.numRows() << std::endl;
    f << " Number of nonzeros = " << pattern.numNonzeros() << std::endl;
    for(unsigned int irow=0;irow < pattern.numRows();irow++){
      f << " Row " << irow << ": ";
      for(int jcol=pattern.rowPointer(irow); jcol < pattern.rowPointer(irow+1);jcol++){
        f << pattern.columnIndex(jcol) << ' ';
      }
      f << '\n';
    }
    f << "......... End of CSR Pattern .........\n";
    return f;
  }

  void CSRMatrixPattern::print(int rankProc) const{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n======== CSR Pattern on Processor #%d  ========\n",rankProc);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Number of rows : %d \n",m_rowptr.size()-1);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Number of nonzeros : %d \n",m_colind.size());
    for(unsigned int irow=0;irow < m_rowptr.size()-1;irow++){
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Row %d: ",irow);
      for(int jcol=m_rowptr[irow]; jcol < m_rowptr[irow+1];jcol++){
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d ",m_colind[jcol]);
      }
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"........ End of CSR Pattern on Processor #%d ........\n",rankProc);

#ifdef PETSC_3_6
    PetscSynchronizedFlush(PETSC_COMM_WORLD,NULL);
#else
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif

  }

}
