#ifndef __CardioXcomp__STIMULATION__
#define __CardioXcomp__STIMULATION__

#include <petsc.h>
#include "mesh.hpp"

#include <vector>
#include "platform.h"

using namespace std;

namespace cardioxcomp{

  class Stimulation{
  protected:

    bool m_isManual; // manual or automatical stimulation
    int m_dim; // dimension of the problem (0D/2D)

    double m_start;
    double m_iapp;
    std::vector<double> m_starts;
    std::vector<double> m_ends;
    std::vector<double> m_iapps;
    double m_repetition; // period of the cycle
    double m_duration; // duration of the stimulation
    
    std::vector<double> m_position; //size 3: circle (xcenter, ycenter, radius), size 4: rectangle (xmin,xmax,ymin,ymax)
    
    std::vector<int> m_nodesInsideArea;

  public:

#ifdef _WIN32
    __declspec(dllexport)
#endif
    Stimulation();

    Stimulation(bool isManual, int dim, const std::vector<double>& starts, const std::vector<double>& ends, const std::vector<double>& iapps);
    Stimulation(bool isManual, int dim, double start, double iapp, double duration, double repetition);

    Stimulation(Mesh* mesh, bool isManual, int dim, const std::vector<double>& starts, const std::vector<double>& ends, const std::vector<double>& iapps, const std::vector<double>& position);
    Stimulation(Mesh* mesh, bool isManual, int dim, double start, double iapp, double duration, double repetition, const std::vector<double>& position);
    
#ifdef _WIN32
    __declspec(dllexport)
#endif	
	~Stimulation();

    
    std::vector<int> nodesInsideArea(Mesh* mesh, const std::vector<double>& focus);
    
    double value_Iapp(int indice, double time);
    bool isInStimulationArea(int indice); // to locate if the node corresponding to the indice is inside the stimulation

  };


}

#endif
