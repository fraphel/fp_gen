#ifndef __CardioXcomp__luafile__
#define __CardioXcomp__luafile__

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <sstream>


extern "C" {

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

}

using namespace std;

namespace cardioxcomp {

  /*
    This class is used to read main variables:
    -Resolution (order, duration, iterations,...)
    -Mesh
    -MEA informations: electrodes,...
    -Stimulation
    -Specific variables for the bidomain equations
  */

  class LuaScript {
  public:
    
    LuaScript(){};
    LuaScript(const string& filename);
    ~LuaScript();

    lua_State* L;////Ajout!

    void printError(const string& variableName, const string& reason);
    std::vector<int> getIntVector(const string& name);
    std::vector<double> getDoubleVector(const string& name);
    std::vector<string> getStringVector(const string& name);
    
    inline void clean();
    
    template<typename T>
    T get(const string& variableName);
    
    bool lua_gettostack(const string& variableName);
    
    // Generic get
    template<typename T>
    T lua_get(const string& variableName){return 0;}
    
    template<typename T>
    T lua_getdefault(){return 0;}
    
  private:
    ///lua_State* L;
    string filename;
    int level;
  };
  
  // Specializations
  template <>
  inline bool LuaScript::lua_get<bool>(const string& variableName){return (bool)lua_toboolean(L, -1);}
  
  template<>
  inline double LuaScript::lua_get<double>(const string& variableName);
  
  template <>
  inline int LuaScript::lua_get<int>(const string& variableName);
  
  template <>
  inline string LuaScript::lua_get<string>(const string& variableName);
  
  template<>
  inline string LuaScript::lua_getdefault<string>(){return "null";}
  
  
  //###########################################################################################################################
  //###########################################################################################################################
  
  /*
    This class is used to read main variables:
     -Resolution (order, duration, iterations,...)
     -Mesh
     -MEA informations: electrodes,...
     -Stimulation
     -Specific variables for the bidomain equations
   */  


  
  
  class LuaFile : public LuaScript{
  public:
    
    LuaFile():LuaScript(){};
    LuaFile(const string& filename);

    void readLua();
    
    // Problem's dimension
    //---------------------------
    
    int dimension;
    
    //---------------------------
    //---------------------------
    
    // Frequency to write solution
    //---------------------------
    
    int freq_writeSolution;
    bool writeEDOs,writeCurrents;
    double timeStartWriting,timeEndWriting;
    bool binaryPostProc;
    
    //---------------------------
    //---------------------------

    // Resolution
    //---------------------------
    
    bool splitting;
    double theta_split;
    
    //CVODE Options
    //---------------------------
    double rtol,atol; //relative and absolute tolerance for CVODE
    int maxOrder,maxNumStep;//CVODE options
    double maxStep;//CVODE options
    int VmOrderBdf;
    
    //Solver Options
    //---------------------------
    int maxIteration;
    double relativeTolerance,absoluteTolerance;
    string solver,preconditioner,preconditionerOption;
    
    //---------------------------
    //---------------------------
    
    // Time informations
    //---------------------------
    
    double dt,tmax;
    int nmax;
    
    //---------------------------
    //---------------------------
    
    // Ionic model
    //---------------------------
    
    string ionicmodel,cells;
    bool spontaneous;
    double Vm_init;
    double Vmin,Vmax;
    bool CVodeSolver;
    //---------------------------
    //---------------------------
    
    // MEA informations
    //---------------------------
    
    string meshfile;
    int nbEl;
    std::vector<double> Ri,Rel,Cel,diametres; 
    double Ue_init;//Field potential value (t=0)
    
    //---------------------------
    //---------------------------
    
    // Boundary conditions
    //---------------------------
    
    std::vector<int> labels;
    std::vector<string> variables;
    std::vector<double> values;
    
    //---------------------------
    //---------------------------
    
    // Stimulation informations
    //---------------------------
    
    bool manual;
    double start,duration,repetition,iapp;
    std::vector<double> starts,ends,iapps,focus;
    
    //---------------------------
    //---------------------------
    
    // Bidomain informations
    //---------------------------
    
    double Am;
    double Cm;
    double sigma_i,sigma_e;
    
    //---------------------------
    //---------------------------

    // Heterogeneity
    //---------------------------

    bool isHeterogeneity;
    bool isAtrialVentricular;
    bool isRescaling;
    std::vector<double> posHeterogeneity;
    std::vector<double> rescalingConstants;

    // CellHeterogeneity
    //---------------------------
    bool cellHeterogeneity;

    //---------------------------
    //---------------------------

    // Compound action
    //---------------------------

    bool isDrug;
    double dose;
    std::vector<string> channel;
    std::vector<double> IC50;
    //bool isActivator;
    //double alpha;
    
    
    
  };

  //_______________________________________________________________________________________________________________________________________________________
  //_______________________________________________________________________________________________________________________________________________________
  
#ifdef _WIN32
  extern "C" __declspec(dllexport) void valeurDim(LuaFile*);
  extern "C" __declspec(dllexport) LuaFile* CreateConfig(int dim, int fwS, bool wODE, bool wCur, bool binPP, bool splitting, double rtol, double atol, int maxOrd, int maxNS, double maxS, int VmOB, double dt, double tmax, int nmax, int maxIt, double reltol, double abstol, char* thesolver, char* thepreconditioner, char* thepreconditionerOption, char* modelionic, double initialVm, double vmin, double vmax, bool CVodeSolver, char* cells, bool spontaneous, char* filemesh, int nbElectrodes, double* deviceResistance, double* electrodeResistance, double* electrodeCapacitance, double* electrodeDiameter, double initialFieldPotential, int* thelabels, int* thevariables, double* thevalues,  bool isManualStim, double* startStim, double* endStim, double* amplitudeStim, int nbStim, double* positionStim, double AmValue, double CmValue, double sigmai, double sigmae, bool isHetero, bool isAV, bool isResc, double* posHetero, int sizeHetero, char* drugName, bool isDrug, double dose, char* channel, double* ic50, int ic50number);
  extern "C" __declspec(dllexport) LuaFile* CreateConfigCpp(int dim, int fwS, bool wODE, bool wCur, bool binPP, bool splitting, double rtol, double atol, int maxOrd, int maxNS, double maxS, int VmOB, double dt, double tmax, int nmax, int maxIt, double reltol, double abstol, string thesolver, string thepreconditioner, string thepreconditionerOption, string modelionic, double initialVm, double vmin, double vmax, bool CVodeSolver, string cells, bool spontaneous, string filemesh, int nbElectrodes, vector<double> deviceResistance, vector<double> electrodeResistance, vector<double> electrodeCapacitance, vector<double> electrodeDiameter, double initialFieldPotential, vector<int> thelabels, vector<string> thevariables, vector<double> thevalues, bool isManualStim, vector<double> startStim, vector<double> endStim, vector<double> amplitudeStim, int nbStim, vector<double> positionStim, double AmValue, double CmValue, double sigmai, double sigmae, bool isHetero, bool isAV, bool isResc, vector<double> posHetero, string name, bool isDrug, double dose, vector<string> channel, vector<double> ic50);
#endif
  
  //_______________________________________________________________________________________________________________________________________________________
  //_______________________________________________________________________________________________________________________________________________________
  
  
  /*
    This class is used to read ionic variables (MV,PCBE,Paci,Courtemanche,Davies):
     -Initials gates values
     -Initials concentrations values
     -Conductances
     -Parameters
  */
  
  class LuaVariables : public LuaScript{
    
  public:
    
    LuaVariables():LuaScript(){};
    LuaVariables(const string& filename);
    
    void readLua(const string& model);
    
    //All initials values (gates+concentrations)
    double hgate,fgate,rgate,sgate,vgate,wgate;
    double jgate,mgate,dgate,fcagate,f1gate,f2gate,qgate,xr1gate,xr2gate,xsgate,xfgate,ggate;
    double Cai,Casr,Nai_v, Nai_a, Nai;
    double mORdgate,jORdgate,hslowgate,hfastgate,hlgate,mlgate;
    

    //O Hara Rudy
    double hCaMKslowgate, jCaMKgate, mLgate, hLgate, hLCaMKgate, agate, ifastgate, islowgate, aCaMKgate;
    double iCaMKfastgate, iCaMKslowgate, ffastgate, fslowgate, fCafastgate, fCaslowgate, jCagate, fCaMKfastgate;
    double fCaCaMKfastgate, ngate, xrfastgate, xrslowgate, xs1gate, xs2gate, xk1gate;
    double JrelNPcurrent, JrelCaMKcurrent, CaMKtrap;
    double Nass, Kss, Cass, Cansr, Cajsr;
    double ead;
    int cellType;

    double oagate,oigate,uagate,uigate,xrgate,ugate;
    double Ki,Caup,Carel;

    //TNNP
    double rprimegate;

    //Davies
    double  dpgate, fca2gate, ydvgate, ydv2gate, zdvgate, Cli, AAgate, rogate, rigate;
    
    //All variables (ionic parameters)
    double f0,beta,eps;
    double V1,V2,Vc,Vs,Vfi;
    double beta1,beta2;
    double tauhp,tauhm,taufp,taufm,taurp,taurm,tausp,tausm;
    double taufCa,taug;
    double tautr,tauu;
    double tauv1m,tauv2m,tauvp,tauw1m,tauw2m,tauwp,taufi,tauo1,tauo2,tauso1,tauso2,taus1,taus2;
    double tausi,tauwinf;
    
    double thetav,thetaw,thetavm,thetao;
    
    double Nao,Cao,Ko,Pkna,Lo,Q,Kmk,KmNa,PNaK,KNaCa,Ksat,KmCa,KmNai,KpCa,Kup,Bufc,Bufsr,Kbufc,Kbufsr;
    double PNaK_v,PNaK_a, KNaCa_v, KNaCa_a;
    double Vmaxup,Vleak,Vsr,alpha,gamma,arel,brel,crel;
    double Vmaxup_v, Vmaxup_a, Vsr_v, Vsr_a, Vc_v, Vc_a;
    double Ef,constf2,V0,Cm;
    double Cm_v, Cm_a, V0_v, V0_a;
    
    double INaKmax,INaCamax,IpCamax,Vi,Vrel,Vup,Caupmax,Iupmax,Cmdnmax,Trpnmax,KmCmdn,KmTrpn,Csqnmax,KmCsqn, Ikurmax;
    double Krel,KmKo;
    
    double u0,uu,kwm,uwm,kso,uso,ks,us,winfstar;
    
    //Conductances
    double gfi,gsi,gto,gso;
    double gNa,gCaL,gKr,gKs,gK1,gf,gbNa,gbCa,gpCa;
    double gNa_v, gNa_a, gto_v, gto_a, gK1_v, gK1_a;
    double gpK;
    double PCa;//ADDED FOR SOBIE ARTICLE
    double gKb, gNa_fast, gNa_late, gNaCa;// for ORd model


    //Assimilation MV
    string drug_gfi, drug_gso, drug_gsi;
    double concentration;
    double configNa, configKCa;
    
    // Davies additions
    double gNaL,gKp,gClb,gto2,Irelmax,Idiffmax,Clo;
    
    std::vector<double> conductances;
    std::vector<double> conductances_v;
    std::vector<double> conductances_a;
    std::vector<double> Imax,V,concentrations,tau,theta,others;
    std::vector<double> V_v,V_a,others_v, others_a;

    std::vector< std::vector<double> > m_InitialCond;
  };

#ifdef _WIN32  
  extern "C" __declspec(dllexport) LuaVariables* CreateVariables(double* initialgates, double* conductances, double* Vparam, double* tauParameters, double* othersParameters, double* concentrationsParameters, char* modelionic);
  extern "C" __declspec(dllexport) LuaVariables* CreateVariablesCpp(string ionicmodel, vector<double> initialgates, vector<double> conductances, vector<double> Vparam, vector<double> tauParameters, vector<double> othersParameters, vector<double> concentrations);
#endif
  
  //###########################################################################################################################
  //###########################################################################################################################
  
  
}

#endif
