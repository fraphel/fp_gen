#ifndef __CardioXcomp__IONICMODEL__
#define __CardioXcomp__IONICMODEL__


#include <vector>
#include <map>
#include <iostream>
#include <math.h>
#include "platform.h"

#include <string>
#include <sstream>

#include <fstream>

#include <petscvec.h>
#include <petscao.h>
#include "dof.hpp"
#include "luafile.hpp"

#include "stimulation.hpp"
#include "heterogeneity.hpp"
#include "mesh.hpp"
//FOR CVODE__________________________________________________________________________________

#include <cvode/cvode.h>          /* prototypes for CVODE fcts. and consts. */
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., and macros */
#include <cvode/cvode_dense.h>    /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DenseMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

//___________________________________________________________________________________________

using namespace std;

namespace cardioxcomp{
  
  static map<string,double> parameters; // global variable containing all ionic parameters
  static Stimulation stimulation;
  static Heterogeneity heterogeneity;
  static bool isHeterogeneity;
  static bool isAtrialVentricular;
  static bool isRescaling;
  
  static double fromSolutionToIonicModel(double,double,double);

  static double m_stockGsi;
  static double m_stockGfi;
  static double m_stockGso;
  static double m_stockGto;

  void particleToVariable(string,double);
  double getVariableValue(string);
  
  //#################################//
  // Data: to communicate with CVODE //
  //#################################//
  
  class UserData {
  public:
    bool m_VmInCVODE;
    bool m_dim0;
    double m_Vm; // useless if Vm is in CVODE else, the value corresponds to the new extrapolate value of Vm
    double m_value_I_ion; // value of I_ion (useless if Vm is in CVODE)
    int m_pos; // useless in 0D
    bool m_atrial; //Is it atrial
    bool m_spontaneous; // for Paci model
    bool m_cellIsRescaled; // for Paci model
    std::vector<double> m_rescalingConstants;
    std::vector<double> m_current_value; // to get the current value if we use a splitting method
    double Vmean;
    int m_numODE;
    double m_vmin,m_vmax, m_Cm;//For phenomenologic model
    double m_heteroCoeff;
    vector<double> m_param;
    bool m_cellHeterogeneity;
    void initializeConductancesPerNode(const int iNode, const double heteroCoeff, const int numCurrents);
    void setCellHeterogeneityFlag(){m_cellHeterogeneity = true;};//(m_cellHeterogeneity = true);
    int m_cellType;
    int m_cellFacto;
  };

  class DrugCompound{
  public:
    bool Drug; // For drug action
    double Dose;
    double hillCoef;  //could be changed to a vector of double to be channel dependent if data is given
    //bool Nablock, Krblock, CaLblock, funnyblock;
    //bool NaKblock, NaCablock, pCablock, bNablock, bCablock, K1block, toblock, Ksblock;
    //double Naic50, Kric50, CaLic50, Ific50;
    //double NaKic50, NaCaic50, pCaic50, bNaic50, bCaic50, K1ic50, toic50, Ksic50;
    std::vector<string> channels;
    std::vector<double> IC50s;
  };

  //#################################//
  
  
  //#######################################################//
  // IonicModel: contains all common variables and methods //
  //#######################################################//
  
  class IonicModel{
  protected:
    double m_time;
    double m_dt;
    int m_numVertex;  // number of vertices in the (global) mesh
    int m_numVertexLocal; // number of vertices in the processor (local size of Vm, ion,...)
    double m_rtol; // for CVODE
    double m_atol; // for CVODE
    int m_maxOrder; // for CVODE
    int m_maxNumStep; // for CVODE
    double m_maxTimeStep; // for CVODE
    int m_numODE; // number of ODEs (including the one for Vm, even it is not solved by CVODE).
    int m_numCurrent; // number of currents in the ionic model
    bool m_VmInCVODE; /* true if Vm is solved by CVODE. In this case the number of equations solved by CVODE is m_numODE.
                         false if Vm is solved elsewhere (typically in the PDE). In this case the number of equations solved by CVODE is m_numODE-1.
                      */
    std::vector<void *> m_cvode_mem; // for CVODE
    std::vector<string> m_listODE; // name of the ODEs
    std::vector<string> m_listCurrent; // name of the ionic currents

    // if many compounds at the same time
    // DrugCompound* Compounds;
    DrugCompound m_compound;
    //===================================================
    
    UserData* m_userData;
    bool m_dim0;
    LuaVariables m_ionicVariables;//user interface
    vector< vector<double> > m_InitialCond;//initial condition if used
    
    realtype m_time_outCVODE; // time returned by CVODE (>= dt) (output time for 0D case)
    
    platformInt m_minOwnershipRange; // global index of the first local vertex
    Vec m_I_ion; // ionic current for the rhs of the PDE
    Vec m_Vm; // 2d case (no splitting)
    std::vector<Vec> m_vecSol; // m_vecSol[i] = petsc vector of the i-th ODE solution (the petsc vector has a size=1 for 0D problems)
    std::vector<Vec> m_vecCurrent; // m_vecCurrent[i] = petsc vector of the i-th current
    std::vector<N_Vector> m_cvodeSol; // m_cvodeSol[i] = cvode vector on the i-th local vertex of the mesh (i=0 for 0D problems)
    std::vector<N_Vector> m_cvodeAbsTol; // m_cvodeAbsTol[i] = cvode absolute tolerance on the i-th local vertex of the mesh (i=0 for 0D problems)
    
  public:
    IonicModel(double t0, double dt, int numVertex, int numVertexLocal, double rtolcvode, double atolcvode, int maxOrdcvode, int maxNStepcvode, double maxTimeStepcvode,int numODE, int numCurrent, bool Vm_in_cvode , bool dim0, string pID="0", bool cellHeterogeneity = false);
    ~IonicModel();  //destructor
    void setMaxTimeStep(double dtmax); // set the maximum time step in CVODE
    void solveODE(double time); // to call CVODE to solve the ODE problem
    void solveODE_RL(double time, double DT, IonicModel *my_model); // to call CVODE to solve the ODE problem, but with the RL integration of the gating var
    void compute_I_ion(); // to compute I_ion and return it in the PDE problem (Vm not in CVODE)
    void actualizeCVODE(double time); // to actualize Vm in CVODE at time t=time (new initial condition on Vm)

    void actualizeCVODE(double time, const vector<double>& ode); // to actualize Vm in CVODE at time t=time (new initial condition on all ODE)

    // Getters ans Setters
    
    bool getVmInCVODE(){return m_VmInCVODE;};
    std::vector<Vec> getCurrent(){return m_vecCurrent;};
    std::vector<Vec> getODE(){return m_vecSol;};
    std::vector<string> getODEname(){return m_listODE;};
    std::vector<string> getCurrentname(){return m_listCurrent;};
    Vec getVm(){return m_Vm;};
    Vec getI_ion(){return m_I_ion;};
    
    int getNVlocal(){return m_numVertexLocal;};
    int getNCurrent(){return m_numCurrent;};
    int getOwnerRange(){return m_minOwnershipRange;};
    std::vector<Vec> getVecSol(){return m_vecSol;};
    void setVecSol(std::vector<Vec> vecSol);
    void setVecCurrent(std::vector<Vec> vecCurr);
    void setI_ion(Vec Iion);
  
    void setVm(Vec Vm);
    double getTime_outCVODE(){return (double)m_time_outCVODE;};
    void setStimulation(const Stimulation& stimulation);
    void setHeterogeneity(const Heterogeneity& heterogeneity);
    void setIsHeterogeneity(bool isHeterogeneity);
    void setIsAtrialVentricular(bool isAtrialVentricular);
    void setIsRescaling(bool isRescaling);
    
    
    void check_variables();
    void check_I_ion();
    void setIonicVariables(const LuaVariables& luavariables);
    
    void readHeterogeneityField(vector<double> &heteroVec, string pID);


    int getInitialCondition(double);
  };
  
  //#######################################################//
  
  //###################//
  // O'Hara Rudy model //
  //###################//
  
  class IonicModel_ORd:
    public IonicModel{

  public:
    IonicModel_ORd(double t0, double dt, int numVertex, int numVertexLocal, double rtolcvode, double atolcvode, int maxOrdcvode, int maxNStepcvode, double maxTimeStepcvode, bool Vm_in_cvode, double Vm0, bool dim0, const LuaVariables& luavar, string pID="0",bool cellHeterogeneity=false);
    static int M_f(realtype t, N_Vector y, N_Vector ydot, void *f_data);   // compute the rhs for the whole system


    static void mgate(double&, double&, double);
    static void hgate(double&, double&, double&, double);
    static void jgate(double&, double&, double);
    static void hCaMKgate(double&, double&, double);
    static void jCaMKgate(double&, double&, double);
    static void mLgate(double&, double&, double);
    static void hLgate(double&, double&, double);
    static void hLCaMKgate(double&, double&, double);
    static void agate(double&, double&, double);
    static void igate(double&, double&, double&, double&, double);
    static void aCaMKgate(double&, double&, double);
    static void iCaMKgate(double&, double&, double&, double&, double);
    static void dgate(double&, double&, double);
    static void fgate(double&, double&, double&, double);
    static void fCagate(double&, double&, double&, double&, double);
    static void jCagate(double&, double&, double);
    static void fCaMKgate(double&, double&, double);
    static void fCaCaMKgate(double&, double&, double);
    static void xrgate(double&, double&, double&, double&, double);
    static void xs1gate(double&, double&, double);
    static void xs2gate(double&, double&, double);
    static void xk1gate(double&, double&, double, double);
  };

  //###################//
  
  
  double Heaviside(double,double);
  
}
#endif
