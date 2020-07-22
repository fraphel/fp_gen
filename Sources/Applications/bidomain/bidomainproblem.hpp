#ifndef __CardioXcomp__bidomainproblem__
#define __CardioXcomp__bidomainproblem__

#ifdef _WIN32
#include "problem.hpp"
#include "bdf.hpp"
#include "platform.h"
#include "ionicmodel.hpp"
#include "matElem.hpp"
#else
#include "../../LibCardioXcomp/problem.hpp"
#include "../../LibCardioXcomp/bdf.hpp"
#include "../../LibCardioXcomp/platform.h"
#include "../../LibCardioXcomp/ionicModel.hpp"
#include "../../LibCardioXcomp/matElem.hpp"
#include "../../LibCardioXcomp/luafile.hpp"
#include "../../LibCardioXcomp/output.hpp"
#endif

#include <petsc.h>
#include "electrodes.hpp"
#include <vector>
#include <petscmat.h>

namespace cardioxcomp {
  class Bidomainproblem:
    public Problem
  {
  public:

    Bidomainproblem();
    Bidomainproblem(std::string meshfile,int numProc,int rankProc, const LuaFile& parameters, const LuaVariables& luavar, const string procId);

    ~Bidomainproblem();

    void initializeProblem(const Vec& uExtrap);
    void initializeElectrodeParameters(int nbEl, const vector<double>& Ri, const vector<double>& Rel, const vector<double>& Cel, const vector<double>& Diameters);
    
    void computeMatVec();
    void dirichletBC();
    
    void setVm_new(); // set the last value of Vm in the ionicModel
    Vec computeODEpart(double time, Vec& V_split, bool actualize);


    void setLuaConfig(const LuaFile& luafile){m_parameters=luafile;};
    void setLuaVar(const LuaVariables& luavar){m_variables=luavar;};
    
    
    void electrodeSimulation(bool actualize, std::vector<double>& Ielsource);
    void electrodeNeumannConditions(const std::vector<double>& Ielsource);
    void computeElectrodesMean();
    
    double computeCaiMean();
    
    void groundBC();
    
    void gatherODE();
    void gatherCurrent();
    void output();
    
    void writeElectrodesSolution();
    void writeElectrodesSolution_headers();
    void writeUemes(const std::vector<double>& Umes);
    void writeUemesAtTheEnd();

    double convertVmForIonicModel(double Vm);
    double convertVmForSolution(double VmScaled);

    int getInitialCondition(const double cell);
    void readHeterogeneityField(vector<double> &heteroVec, string pID);
      
    LuaFile m_parameters;
    LuaVariables m_variables;


    Vec m_ion; // = I_ion + Electrode (if no splitting), = Electrode (if splitting)
    
    //Ionic model
    IonicModel* m_ionicModel;

    string m_ionicmodel;
    bool m_atrial,m_spontaneous;

    double m_rtol,m_atol;
    int m_maxOrder,m_maxNumStep;
    double m_maxStep;
    double m_sigmai,m_sigmae;
    double m_Am,m_Cm;//problem (rhs)

    Vec m_vecVmForIonicModel;
    double m_vmin, m_vmax, m_dv;

    //########//
    // Output //
    //########//

    bool m_writeODE,m_writeCurrent;
    string m_procId;

    //########//
    //########//

    int getNbEl(){ return m_nbEl;};

    //##########//
    // For Iapp //
    //##########//

    Stimulation m_stim;
    
    //##########//
    //##########//

    //###############//
    // Heterogeneity //
    //###############//

    Heterogeneity m_heterogeneity;
    bool m_isHeterogeneity;
    bool m_isAtrialVentricular;
    bool m_isRescaling;
    
    
    //####################//
    // Cell Heterogeneity //
    //####################//
    bool m_cellHeterogeneity;


    //############//
    // Electrodes //
    //############//
    
    int m_nbEl;
    std::vector<double> m_Uemean,m_Vmmean;//Moyenne aux electrodes
    double m_CaimeanWell;
    vector< double > m_CaiVec;
    vector< vector<double> > m_UmesVec, m_VmVec;
    std::vector<double> m_measure;
    std::vector<int> m_labels;
    Electrodes m_electrodes;
    std::vector<double> m_Ri,m_Rel,m_Cel,m_diametres;

    //############//
    //############//
    
    bool m_isCVOdeSol;
    
    //####################//
    // Boundary Condition //
    //####################//
    
    std::vector<int> m_labelsBC;
    std::vector<string> m_variablesBC;
    std::vector<double> m_valuesBC;

    //####################//
    //####################//
    
    
    //##########//
    // Postproc //
    //##########//
    
    std::vector<string> m_list_ODE,m_list_Current;
    std::vector<Vec> m_vecODE,m_vecCurrent;
    std::vector<Vec> m_vecSeqODE,m_vecSeqCurrent;
    
    //##########//
    //##########//
    
    //void check_variables(std::vector<Vec> vecODE);
    
  };
}

#endif

