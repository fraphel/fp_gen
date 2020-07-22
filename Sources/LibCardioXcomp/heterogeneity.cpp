#include "stdafx.h"
#include "heterogeneity.hpp"



namespace cardioxcomp{

#ifdef _WIN32
  __declspec(dllexport)
#endif
  Heterogeneity::Heterogeneity(){

  }

  Heterogeneity::Heterogeneity(Mesh* mesh, const std::vector<double>& area){
    m_nodesInsideHetero = nodesInsideHetero(mesh,area);
  }

#ifdef _WIN32
  __declspec(dllexport)
#endif
  Heterogeneity::~Heterogeneity(){
  }



  std::vector< std::vector<int> > Heterogeneity::nodesInsideHetero(Mesh* mesh, const std::vector<double>& area){

    std::vector< std::vector<int> > hetero;

    int nb_hetero = 1;
    hetero.resize(nb_hetero);

    for(unsigned int ih=0 ; ih<hetero.size() ; ih++){
      
      for(int indice=0 ; indice<mesh->nv ; indice++){
        
        if(mesh->vertices[indice].x >= area[0] && mesh->vertices[indice].x <= area[1] && mesh->vertices[indice].y >= area[2] && mesh->vertices[indice].y <= area[3]){
          hetero[ih].push_back(indice);
        }
        
      }
      
    }
    
    return hetero;
  }


  bool Heterogeneity::isInHetero(int indice){
    
    bool isInside = false;
    
    for(unsigned int inode=0 ; inode<m_nodesInsideHetero[0].size() ; inode++){
      if(indice == m_nodesInsideHetero[0][inode]){
        isInside = true;
        break;
      }
    }
    
    return isInside;
  }
  


}
