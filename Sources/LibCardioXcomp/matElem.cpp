#include "stdafx.h"

#include "matElem.hpp"


namespace cardioxcomp {
  
  void MatElem::init(int numComp){
    m_numComp=numComp;
    m_numComp2=numComp*numComp;
    m_size=9*m_numComp2;
    val = new double[m_size];
  }

  void MatElem::zero(){
    for(int i=0;i<m_size;i++) val[i] = 0.;
  }

}
