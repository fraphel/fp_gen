#include "stdafx.h"

#include "stimulation.hpp"



namespace cardioxcomp{

#ifdef _WIN32
  __declspec(dllexport)
#endif
  Stimulation::Stimulation(){
  }
  
#ifdef _WIN32
  __declspec(dllexport)
#endif
  Stimulation::~Stimulation(){
  }
  
  
  //############//
  // Manual, 0D //
  //############//
  
  Stimulation::Stimulation(bool isManual, int dim, const std::vector<double>& starts, const std::vector<double>& ends, const std::vector<double>& iapps){
    
    m_isManual = isManual;
    m_dim = dim;
    m_starts = starts;
    m_ends = ends;
    m_iapps = iapps;
    
  }
  
  
  //############//
  // Manual, 2D //
  //############//
  
  Stimulation::Stimulation(Mesh* mesh, bool isManual, int dim, const std::vector<double>& starts, const std::vector<double>& ends, const std::vector<double>& iapps, const std::vector<double>& position){
    
    m_isManual = isManual;
    m_dim = dim;
    m_starts = starts;
    m_ends = ends;
    m_iapps = iapps;
    
    m_position = position;
    
    m_nodesInsideArea = nodesInsideArea(mesh,position);
    
  }
  
  
  //###############//
  // Automatic, 0D //
  //###############//
  
  Stimulation::Stimulation(bool isManual, int dim, double start, double iapp, double duration, double repetition){
    
    m_isManual = isManual;
    m_dim = dim;
    m_start = start;
    m_iapp = iapp;
    
    m_duration = duration;
    m_repetition = repetition;
    
  }
  
  
  //###############//
  // Automatic, 2D //
  //###############//
  
  Stimulation::Stimulation(Mesh* mesh, bool isManual, int dim, double start, double iapp, double duration, double repetition, const std::vector<double>& position){
    
    m_isManual = isManual;
    m_dim = dim;
    m_start = start;
    m_iapp = iapp;
    
    m_duration = duration;
    m_repetition = repetition;
    
    m_position = position;
    
    m_nodesInsideArea = nodesInsideArea(mesh,position);
    
  }
  
  
  double Stimulation::value_Iapp(int indice, double time){
    
    double iapp = 0.;
    
    if(m_isManual){
      
      for(unsigned int istim=0 ; istim<m_starts.size() ; istim++){
        
        if( (m_starts.size() != m_ends.size()) || (m_ends.size() != m_iapps.size()) ){
          std::cout << "ERROR in stimulation parameters: start, end and iapp must have the same length" << std::endl;
          exit(0);
        }
        
        if(time>=m_starts[istim] && time<=m_ends[istim]){
          if(m_dim==0){
            iapp = m_iapps[istim];
          }
          else{
            if(Stimulation::isInStimulationArea(indice)){
              iapp = m_iapps[istim];
            }
          }
        }
        else if(time<=m_starts[istim] && time>=m_ends[istim]){
          std::cout << "WARNING: Check the time stimulation"<< std::endl;
        }
        
      }
    }
    else{
      
      const int numStim = time/m_repetition;
//      std::cout << "time = " << time << ", numStim = " << numStim << std::endl;
      if( (time>= m_start+numStim*m_repetition) && (time<= m_start+numStim*m_repetition+m_duration) ){
        if(m_dim==0){
          iapp = m_iapp;
                    //          std::cout << "       rep = " << m_repetition << " iapp = " << iapp << std::endl;
        }
        else{
          double t0 = m_start+numStim*m_repetition+.5*m_duration;
          if(Stimulation::isInStimulationArea(indice)){
            double sigma = m_duration/6.;
            iapp = m_iapp*exp(-.5*pow((time-t0),2)/pow(sigma,2));
            //iapp = m_iapp;
          }
        }
      }
      
    }
    
    return iapp;
  }
  
  
  bool Stimulation::isInStimulationArea(int indice){
    
    bool isInside = false;
    
    for(unsigned int inode=0 ; inode<m_nodesInsideArea.size() ; inode++){
      if(indice == m_nodesInsideArea[inode]){
        isInside = true;
        break;
      }
    }
    
    return isInside;
  }
  
  
  std::vector<int> Stimulation::nodesInsideArea(Mesh* mesh, const std::vector<double>& focus){
    
    std::vector<int> indices;
    
    for(int inode=0 ; inode<mesh->nv ; inode++){
      
      int indice = inode;
      
      if(focus.size() == 3){
        
        const double r2 = pow(focus[2],2);
        const double x1 = pow(mesh->vertices[indice].x-focus[0],2);
        const double y1 = pow(mesh->vertices[indice].y-focus[1],2);
        const double c0 = x1+y1;
        
        if(c0<r2){
          indices.push_back(indice);
        }
        
      }
      else if(focus.size() == 4){
          
          if(mesh->vertices[indice].x >= focus[0] && mesh->vertices[indice].x <= focus[1] && mesh->vertices[indice].y >= focus[2] && mesh->vertices[indice].y <= focus[3]){
            indices.push_back(indice);
          }
        
        }
        else{
          std::cout << "Error in the configuration lua file (Stimulation->Position): size 3: (xcenter,ycenter,radius) for a circle stimulation, size 4: (xmin,xmax,ymin,ymax) for a rectangular stimulation" << std::endl;
          exit(0);
        }
        
      }
      
      return indices;
    }
    
  }
