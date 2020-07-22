#include "stdafx.h"

#include "vecElem.hpp"

namespace cardioxcomp {
  
  void VecElem::init(int numComp){
    m_numComp = numComp;
    m_size = 3*numComp;
    val = new double[m_size];
  }

  void VecElem::zero(){
    for(int i=0;i<m_size;i++) val[i] = 0.;
  }

}
