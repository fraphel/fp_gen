#ifndef __CardioXcomp__HETEROGENEITY__
#define __CardioXcomp__HETEROGENEITY__

#include <petsc.h>
#include "mesh.hpp"
#include <vector>

using namespace std;

namespace cardioxcomp{

  class Heterogeneity{

  public:

#ifdef _WIN32
    __declspec(dllexport)
#endif
    Heterogeneity();

    Heterogeneity(Mesh* mesh, const std::vector<double>& area);

#ifdef _WIN32
    __declspec(dllexport)
#endif

    ~Heterogeneity();

    std::vector< std::vector<int> > nodesInsideHetero(Mesh* mesh, const std::vector<double>& area);
    bool isInHetero(int indice);


  protected:
    std::vector< std::vector<int> > m_nodesInsideHetero;

  };


}



#endif
