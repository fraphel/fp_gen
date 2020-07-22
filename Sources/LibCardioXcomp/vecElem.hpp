#ifndef __CardioXcomp__vecelem__
#define __CardioXcomp__vecelem__

namespace cardioxcomp {
  class VecElem{
    int m_numComp;
    int m_size;
  public:
    VecElem(){}
    void init(int numComp);
    double* val;
    void zero();
    inline int index(int icomp,int i){return i*m_numComp + icomp;}
  };
}

#endif
