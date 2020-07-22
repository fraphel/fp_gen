#ifndef __CardioXcomp__electrodes__
#define __CardioXcomp__electrodes__

#include <vector>
#ifdef _WIN32
#include "bdf.hpp"
#include "platform.h"
#else
#include "../../LibCardioXcomp/bdf.hpp"
#include "../../LibCardioXcomp/platform.h"
#endif

using namespace std;


namespace cardioxcomp{

  class Electrodes{

  public:
    
    Electrodes();
    ~Electrodes();
    
    void initializeElectrodes(int nE, double dt, const vector<double>& C, const vector<double>& Ri, const vector<double>& Rel, const vector<double>& diameters);
    void computeAreas();
    void computeTau();
    void computeRtilde();
    void forwardElectrodes();
    void computedU();
    void set_Umean(const vector<double>& umean);
    void computeUmes();
    void getUmes(vector<double>&);
    void computeSourceIel(vector<double>&);

    void computeIel();

    int m_nbElectrodes;
    double m_dt;

    vector<double> m_C;// capacitance
    vector<double> m_Rel, m_Ri;// resistances
    vector<double> m_diameters;

    vector<double> m_areas;
    vector<double> m_tauk;//C*(Ri+Rel)
    vector<double> m_rtilde;//(Ri+Rel)

    vector<double> m_Umean,m_dU;
    vector<double> m_Umes;

    Bdf m_bdfIel;
    Bdf m_bdfUe;


    Vec m_Iel_Extrap,m_Ue_Extrap;
    Vec m_Iel,m_Ue;
    Vec dU;

  };


}



#endif
