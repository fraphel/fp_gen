#ifndef __CardioXcomp__matelem__
#define __CardioXcomp__matelem__


namespace cardioxcomp {
  class MatElem{
    int m_numComp;
    int m_numComp2;
    int m_size;
  public:
    MatElem(){}
    void init(int numComp);
    void zero();
    double* val;
    inline int index(int icomp,int jcomp,int i,int j){return 3*m_numComp2*i + 3*m_numComp*icomp + j*m_numComp + jcomp;}
  };
}

#endif
