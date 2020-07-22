#include "stdafx.h"

#include "ionicModel.hpp"
#include <iomanip>

using namespace std;
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


#define Rm 8.314472 // [J/mol/K] Gas constant
#define Tcell 310.0 // [K] Temperature

#define Frdy 96485.3415 // Paci, O Hara Rudy, TNNP
#define frdy 96.4867 // Courtemanche (not same unit)
#define RmOHaraRudy 8314.0 // O Hara Rudy
#define RmTNNP 8314.0 // TNNP
#define CmOHaraRudy  1.0
#define CmTNNP 0.185

namespace cardioxcomp{




  void UserData::initializeConductancesPerNode(const int iNode, const double heteroCoeff, const int numCurrents){
    
    /* For OHara Rudy model */
    if (numCurrents == 40){
      if (iNode==0){
        cout << "Ohara-Rudy ionic model. Warning: taking number of currents for heterogeneity mode selection... This must be changed." << endl;
      }
      if (heteroCoeff<1./3.){
        m_cellType = 0;
      }
      else if (heteroCoeff>2./3.){
        m_cellType = 2;
      }
      else{
        m_cellType = 1;
      }
    }
    
  }


  void IonicModel::readHeterogeneityField(vector<double> &heteroVec, string pID){
    // Read heterogeneous conductances
    stringstream fileName;
    fileName << "./hetero/heterogeneity_p" << pID << ".txt"; 
    ifstream inputFile(fileName.str().c_str());
    if (!inputFile){
      cout << "Error. File " << fileName.str() << " does not exist." <<endl;
      exit(1);
    }
    heteroVec.clear();
    double val;
    while (true){
      inputFile >> val;
      heteroVec.push_back(val);
      if (inputFile.eof()) break;
    }
    inputFile.close();
  }




#pragma region Others methods

  void particleToVariable(string name, double value)
  {
    parameters[name] = value;
  }
  
  double getVariableValue(string name)
  {

    double value = parameters[name];
    return value;
  }

  void IonicModel::setIonicVariables(const LuaVariables& luavariables){
	  this->m_ionicVariables = luavariables;
  }

  double Heaviside(double V,double Vsubstr) {
    
    double H = 0.;
    const double substr = V-Vsubstr;
    
    if(substr>=0.) {
      H = 1.;
    }
    return H;
    //return .5 + .5*tanh(1.0e2*(V-Vsubstr));
  }


  
  int IonicModel::getInitialCondition(double cell){
    if (cell>0.3 and cell<0.7){
      return 0;
    }
    else{
      return 1;
    }
  }

  
#pragma endregion

#pragma region IonicModel constructor
  IonicModel::IonicModel(double t0, double dt, int numVertex, int numVertexLocal,double rtolcvode, double atolcvode, int maxOrdcvode, int maxNStepcvode, double maxTimeStepcvode, int numODE, int numCurrent, bool Vm_in_cvode, bool dim0, string pID, bool cellHeterogeneity):
    m_time(t0),
    m_dt(dt),
    m_numVertex(numVertex),
    m_numVertexLocal(numVertexLocal),
    m_rtol(rtolcvode),
    m_atol(atolcvode),
    m_maxOrder(maxOrdcvode),
    m_maxNumStep(maxNStepcvode),
    m_maxTimeStep(maxTimeStepcvode),
    m_numODE(numODE),
    m_numCurrent(numCurrent),
    m_VmInCVODE(Vm_in_cvode),
    m_dim0(dim0)
  {
    	
    m_cvode_mem.resize(m_numVertexLocal);
    m_userData = new UserData[m_numVertexLocal];

    // vector for the ionic current in each node
    VecCreate(PETSC_COMM_WORLD,&m_I_ion);
    VecSetSizes(m_I_ion,m_numVertexLocal,m_numVertex);
    VecSetFromOptions(m_I_ion);
    VecZeroEntries(m_I_ion);

    // vector for the transmembrane potential in each node
    VecCreate(PETSC_COMM_WORLD,&m_Vm);
    VecSetSizes(m_Vm,m_numVertexLocal,m_numVertex);
    VecSetFromOptions(m_Vm);
    VecZeroEntries(m_Vm);
    
    platformInt dummy;
    VecGetOwnershipRange(m_I_ion,&m_minOwnershipRange,&dummy);
    
    if(m_VmInCVODE){
      m_vecSol.resize(m_numODE);
    }
    else{
      m_vecSol.resize(m_numODE-1); // since Vm is not solved in CVODE
    }
    
    // create the vectors for the solutions
    for(unsigned int iODE=0 ; iODE<m_vecSol.size() ; iODE++){
      VecCreate(PETSC_COMM_WORLD,&m_vecSol[iODE]);
      VecSetSizes(m_vecSol[iODE],m_numVertexLocal,m_numVertex);
      VecSetFromOptions(m_vecSol[iODE]);
    }
    // create vectors for the currents
    m_vecCurrent.resize(m_numCurrent);
    for(unsigned int i=0 ; i<m_vecCurrent.size() ; i++){
      VecCreate(PETSC_COMM_WORLD,&m_vecCurrent[i]);
      VecSetSizes(m_vecCurrent[i],m_numVertexLocal,m_numVertex);
      VecSetFromOptions(m_vecCurrent[i]);
      VecZeroEntries(m_vecCurrent[i]);  // initialize currents to zero
    }
    
    
    m_cvodeSol.resize(m_numVertexLocal);
    m_cvodeAbsTol.resize(m_numVertexLocal);

    vector<double> heteroVec;
    if (cellHeterogeneity){
      cout << "THERE IS CELL HETEROGENEITY IN THIS WELL." << endl << "Reading the heterogeneity file..."; 
      readHeterogeneityField(heteroVec, pID);
      cout << "Done." << endl;
    }
    
    for(int i=0 ; i<m_numVertexLocal; i++){
      
      m_cvode_mem[i] = NULL;
      
      m_cvodeSol[i] = N_VNew_Serial(m_vecSol.size());   // vectors (one for each node) with the variables unknown 
      m_cvodeAbsTol[i] = N_VNew_Serial(m_vecSol.size());  // vectors (one for each node) with the tolerances 
      
      m_cvode_mem[i] = CVodeCreate(CV_BDF,CV_NEWTON);  // initiate the solver (stiff problem)
      
      m_userData[i].m_cellHeterogeneity = true;//false;
      if (cellHeterogeneity){
        //m_userData[i].initializeConductancesPerNode(i, heteroVec[i], m_numCurrent);
        //m_userData[i].setCellHeterogeneityFlag();
      }
    }
    
  }
  
  void IonicModel::setMaxTimeStep(double dtmax){
    for(int i=0 ; i<m_numVertexLocal; i++){
      int flag = CVodeSetMaxStep(m_cvode_mem[i],dtmax);
    }
  }

#pragma endregion

#pragma region IonicModel destructor
  IonicModel::~IonicModel(){

    for(unsigned int iODE=0 ; iODE<m_vecSol.size() ; iODE++){
      VecDestroy(&m_vecSol[iODE]);
    }

    for(unsigned int ic=0 ; ic<m_vecCurrent.size() ; ic++){
      VecDestroy(&m_vecCurrent[ic]);
    }

    VecDestroy(&m_I_ion);
    VecDestroy(&m_Vm);

    delete[] m_userData;

  }
#pragma endregion

#pragma region solveODE

  void IonicModel::solveODE(double time){
    
    for(int i=0 ; i<m_numVertexLocal ; i++){  // loop over the nodes
      
      realtype tout;
      platformInt posglob = m_minOwnershipRange+i;
      
      double vm;
      VecGetValues(m_Vm,1,&posglob,&vm);
      m_userData[i].m_Vm = vm;

      int flag = 0;
      // integration of the ODE system:
      /*
       'time' is the next time at which a computed solution is desired (input)
       m_cvodeSol[i] is the computed solution vector.
       tout is the time reached by the solver (output)
       */
      flag = CVode(m_cvode_mem[i], time, m_cvodeSol[i], &tout, CV_NORMAL);
      //std::cout << "### time = " << time << " tout = " << tout << std::endl;
      if(flag<0){
        std::cout << "ionicModel.cpp: CVOde failure" << endl;
        exit(1);
      }
      flag = CVodeGetCurrentTime(m_cvode_mem[i], &m_time_outCVODE);
      //std::cout << " m_time_outCVODE = " << m_time_outCVODE << std::endl;
      
      for(size_t iODE=0 ; iODE<m_vecSol.size() ; iODE++) {  // saving values of the solution inside the vector
        const double valuem_solODE = Ith(m_cvodeSol[i],iODE);
        VecSetValues(m_vecSol[iODE],1,&posglob,&valuem_solODE, INSERT_VALUES);
      }
      
      for(unsigned int ic=0 ; ic<m_vecCurrent.size() ; ic++){  // saving values of the currents inside the vector
        VecSetValues(m_vecCurrent[ic],1,&posglob,&m_userData[i].m_current_value[ic], INSERT_VALUES);
      }

      
    }
    
    // assemble vectors
    for(size_t iODE=0 ; iODE<m_vecSol.size() ; iODE++) {
      VecAssemblyBegin(m_vecSol[iODE]);
      VecAssemblyEnd(m_vecSol[iODE]);
    }

    for(unsigned int ic=0 ; ic<m_vecCurrent.size() ; ic++){
      VecAssemblyBegin(m_vecCurrent[ic]);
      VecAssemblyEnd(m_vecCurrent[ic]);
    }
    
    
  }

  void IonicModel::check_I_ion(){
  
    for(int i=0 ; i<m_numVertexLocal ; i++){  // loop over the nodes
      
      platformInt posglob = m_minOwnershipRange+i;
      
      // saving values of the solution inside the vector
      double valuem_solODE;
      VecGetValues(m_I_ion,1,&posglob,&valuem_solODE);
      if(valuem_solODE > 15000.){
        std::cout<<"PROBLEM!"<<std::endl;
        std::cout<<"Iion "<<" too big"<<std::endl;}
#ifdef _WIN32 //C++ norm not the same
      if(_isnan(valuem_solODE)){
#else
      if(isnan(valuem_solODE)){
#endif
        std::cout<<"PROBLEM!"<<std::endl;
        std::cout<<"Variable "<<" NaN"<<std::endl;}
#ifdef _WIN32 //C++ norm not the same
      if(!_finite(valuem_solODE)){
#else
      if(!isfinite(valuem_solODE)){
#endif
        std::cout<<"PROBLEM!"<<std::endl;
        std::cout<<"Variable "<<" infinite"<<std::endl;}
    }

  }
  

  void IonicModel::check_variables(){
  
    for(int i=0 ; i<m_numVertexLocal ; i++){  // loop over the nodes
      
      platformInt posglob = m_minOwnershipRange+i;
      
      for(size_t iODE=0 ; iODE<m_vecSol.size() ; iODE++) {  // saving values of the solution inside the vector
        double valuem_solODE;
        VecGetValues(m_vecSol[iODE],1,&posglob,&valuem_solODE);
        if(valuem_solODE > 15.){
          std::cout<<"PROBLEM!"<<std::endl;
          std::cout<<"Variable "<<iODE<<" too big"<<std::endl;}
#ifdef _WIN32
        if(_isnan(valuem_solODE)){
#else
        if(isnan(valuem_solODE)){
#endif
          std::cout<<"PROBLEM!"<<std::endl;
          std::cout<<"Variable "<<iODE<<" NaN"<<std::endl;}
#ifdef _WIN32
		if(!_finite(valuem_solODE)){
#else
        if(!isfinite(valuem_solODE)){
#endif
          std::cout<<"PROBLEM!"<<std::endl;
          std::cout<<"Variable "<<iODE<<" infinite"<<std::endl;}
        
      }
    }
 
  }
  
#pragma endregion

#pragma region compute I_ion

  void IonicModel::compute_I_ion(){
    
    VecZeroEntries(m_I_ion);
    
    for(int i=0 ; i<m_numVertexLocal ; i++){
      
      platformInt posglob = m_minOwnershipRange+i;      
      const double value_I_ion = m_userData[i].m_value_I_ion;
      VecSetValues(m_I_ion,1,&posglob,&value_I_ion, INSERT_VALUES);
      
      for(unsigned int ic=0 ; ic<m_userData[i].m_current_value.size() ; ic++){
        VecSetValues(m_vecCurrent[ic],1,&posglob,&m_userData[i].m_current_value[ic], INSERT_VALUES);
      }
    }
    
    
    VecAssemblyBegin(m_I_ion);
    VecAssemblyEnd(m_I_ion);
    
    for(unsigned int ic=0 ; ic<m_vecCurrent.size() ; ic++){
      VecAssemblyBegin(m_vecCurrent[ic]);
      VecAssemblyEnd(m_vecCurrent[ic]);
    }
    
  }
#pragma endregion

#pragma region Actualize CVODE (for splitting)
  void IonicModel::actualizeCVODE(double time){
    
    for(int i=0 ; i<m_numVertexLocal ; i++){
      
      double value;
      platformInt idofloc = i+m_minOwnershipRange;
      VecGetValues(m_Vm,1,&idofloc,&value);
      
      Ith(m_cvodeSol[i],m_numODE-1) = value;      
      CVodeReInit(m_cvode_mem[i],time,m_cvodeSol[i]);//Reinitialize CVODE (new initial condition)
      
    }
    
  }
#pragma endregion

#pragma region Same but for a vector
  void IonicModel::actualizeCVODE(double time, const vector<double>& ode){
    
    for(int i=0 ; i<m_numVertexLocal ; i++){
      
      for(int ie=0 ; ie<m_numODE ; ie++){
        Ith(m_cvodeSol[i],ie) = ode[ie];
      }

      CVodeReInit(m_cvode_mem[i],time,m_cvodeSol[i]);//Reinitialize CVODE (new initial condition)
      
    }
    
  }
#pragma endregion

#pragma region Set new Vm value
  void IonicModel::setVm(Vec Vm){
    VecCopy(Vm,m_Vm); // m_Vm = Vm
  }
#pragma endregion
  
  void IonicModel::setVecSol(std::vector<Vec> vecSol){
    for(size_t iODE=0 ; iODE<vecSol.size() ; iODE++) {  
      VecCopy(vecSol[iODE],m_vecSol[iODE]);  // m_vecSol[i] = vecSol[i]
    }
  }
  
  void IonicModel::setVecCurrent(std::vector<Vec> vecCurr){
    for(size_t iODE=0 ; iODE<vecCurr.size() ; iODE++) {  
      VecCopy(vecCurr[iODE],m_vecCurrent[iODE]);  // m_vecCurrent[i] = vecCurr[i]
    }
  }
  void IonicModel::setI_ion(Vec Iion){
    VecCopy(Iion,m_I_ion); // m_I_ion = Iion
  }
  
#pragma region Stimulation & Heterogeneity
  void IonicModel::setStimulation(const Stimulation& stim){
    stimulation = stim;
  }

  void IonicModel::setHeterogeneity(const Heterogeneity& hetero){
    heterogeneity = hetero;
    // Also modify m_userData accordingly
    for(int i=0 ; i<m_numVertexLocal ; i++)
      {
    	m_userData[i].m_atrial = heterogeneity.isInHetero(m_minOwnershipRange+i);
      }
  }

  void IonicModel::setIsHeterogeneity(bool isHetero){
    isHeterogeneity = isHetero;
  }
#pragma endregion
  
  void IonicModel::setIsAtrialVentricular(bool isHAV){
    isAtrialVentricular = isHAV;
  }
  
  void IonicModel::setIsRescaling(bool isHR){
    isRescaling = isHR;
  }
  
  static double fromSolutionToIonicModel(double vm, double vmin, double vmax){

    double vmIm = 0.;

    vmIm = (vm-vmin)/(vmax-vmin);
    return vmIm;

  }
  
  //__________________________________________________________________________________________________________________________________________________________________________
  //__________________________________________________________________________________________________________________________________________________________________________
  
  // O'Hara Rudy model
  
  //__________________________________________________________________________________________________________________________________________________________________________
  //__________________________________________________________________________________________________________________________________________________________________________
  
#pragma region ORd model
  IonicModel_ORd::IonicModel_ORd(double t0, double dt, int numVertex, int numVertexLocal, double rtolcvode, double atolcvode, int maxOrdcvode, int maxNStepcvode, double maxTimeStepcvode, bool Vm_in_cvode, double Vm0, bool dim0, const LuaVariables& luavar, string pID,bool cellHeterogeneity):
    IonicModel(t0,dt,numVertex,numVertexLocal,rtolcvode,atolcvode,maxOrdcvode,maxNStepcvode,maxTimeStepcvode,41,40,Vm_in_cvode,dim0,pID,cellHeterogeneity)
  {
    std::cout << "CONSTRUCTEUR IonicModel_ORd (O Hara Rudy) \n";
    
    m_listODE.push_back("m_gate");
    m_listODE.push_back("hfast_gate");
    m_listODE.push_back("hslow_gate");
    m_listODE.push_back("j_gate");
    m_listODE.push_back("hCaMK_gate");
    m_listODE.push_back("jCaMK_gate");
    m_listODE.push_back("mL_gate");
    m_listODE.push_back("hL_gate");
    m_listODE.push_back("hLCaMK_gate");
    m_listODE.push_back("a_gate");
    m_listODE.push_back("ifast_gate");
    m_listODE.push_back("islow_gate");
    m_listODE.push_back("aCaMK_gate");
    m_listODE.push_back("iCaMKfast_gate");
    m_listODE.push_back("iCaMKslow_gate");
    m_listODE.push_back("d_gate");
    m_listODE.push_back("ffast_gate");
    m_listODE.push_back("fslow_gate");
    m_listODE.push_back("fCafast_gate");
    m_listODE.push_back("fCaslow_gate");
    m_listODE.push_back("jCa_gate");
    m_listODE.push_back("fCaMKfast_gate");
    m_listODE.push_back("fCaCaMKfast_gate");
    m_listODE.push_back("n_gate");
    m_listODE.push_back("xrfast_gate");
    m_listODE.push_back("xrslow_gate");
    m_listODE.push_back("xs1_gate");
    m_listODE.push_back("xs2_gate");
    m_listODE.push_back("xk1_gate");
    m_listODE.push_back("CaMKtrap");
    m_listODE.push_back("JrelNP_cur");
    m_listODE.push_back("JrelCaMK_cur");
    m_listODE.push_back("Nai");
    m_listODE.push_back("Nass");
    m_listODE.push_back("Ki");
    m_listODE.push_back("Kss");
    m_listODE.push_back("Cai");
    m_listODE.push_back("Cass");
    m_listODE.push_back("Cansr");
    m_listODE.push_back("Cajsr");
    if(Vm_in_cvode) m_listODE.push_back("Vm");
    
    m_listCurrent.push_back("INa");
    m_listCurrent.push_back("Ito");
    m_listCurrent.push_back("ICaL");
    m_listCurrent.push_back("ICaNa");
    m_listCurrent.push_back("ICaK");
    m_listCurrent.push_back("IKr");
    m_listCurrent.push_back("IKs");
    m_listCurrent.push_back("IK1");
    m_listCurrent.push_back("INaCa_i");
    m_listCurrent.push_back("INaCa_ss");
    m_listCurrent.push_back("INaCa");
    m_listCurrent.push_back("INaK");
    m_listCurrent.push_back("INab");
    m_listCurrent.push_back("ICab");
    m_listCurrent.push_back("IKb");
    m_listCurrent.push_back("IpCa");//16
    
    //and fluxes....
    m_listCurrent.push_back("JdiffNa");
    m_listCurrent.push_back("JdiffK");
    m_listCurrent.push_back("JdiffCa");
    m_listCurrent.push_back("JrelNP");
    m_listCurrent.push_back("JrelCaMK");
    m_listCurrent.push_back("Jrel");
    m_listCurrent.push_back("JupNP");
    m_listCurrent.push_back("JupCaMK");
    m_listCurrent.push_back("Jup");
    m_listCurrent.push_back("Jtr");
    
    
    m_ionicVariables = luavar;
    parameters["m0"] = m_ionicVariables.mgate;
    parameters["hfast0"] = m_ionicVariables.hfastgate;
    parameters["hslow0"] = m_ionicVariables.hslowgate;
    parameters["j0"] = m_ionicVariables.jgate;
    parameters["hCaMKslow0"] = m_ionicVariables.hCaMKslowgate;
    parameters["jCaMK0"] = m_ionicVariables.jCaMKgate;
    parameters["mL0"] = m_ionicVariables.mLgate;
    parameters["hL0"] = m_ionicVariables.hLgate;
    parameters["hLCaMK0"] = m_ionicVariables.hLCaMKgate;
    parameters["a0"] = m_ionicVariables.agate;
    parameters["ifast0"] = m_ionicVariables.ifastgate;
    parameters["islow0"] = m_ionicVariables.islowgate;
    parameters["aCaMK0"] = m_ionicVariables.aCaMKgate;
    parameters["iCaMKfast0"] = m_ionicVariables.iCaMKfastgate;
    parameters["iCaMKslow0"] = m_ionicVariables.iCaMKslowgate;
    parameters["d0"] = m_ionicVariables.dgate;
    parameters["ffast0"] = m_ionicVariables.ffastgate;
    parameters["fslow0"] = m_ionicVariables.fslowgate;
    parameters["fCafast0"] = m_ionicVariables.fCafastgate;
    parameters["fCaslow0"] = m_ionicVariables.fCaslowgate;
    parameters["jCa0"] = m_ionicVariables.jCagate;
    parameters["fCaMKfast0"] = m_ionicVariables.fCaMKfastgate;
    parameters["fCaCaMKfast0"] = m_ionicVariables.fCaCaMKfastgate;
    parameters["n0"] = m_ionicVariables.ngate;
    parameters["xrfast0"] = m_ionicVariables.xrfastgate;
    parameters["xrslow0"] = m_ionicVariables.xrslowgate;
    parameters["xs10"] = m_ionicVariables.xs1gate;
    parameters["xs20"] = m_ionicVariables.xs2gate;
    parameters["xk10"] = m_ionicVariables.xk1gate;
    parameters["CaMKtrap0"] = m_ionicVariables.CaMKtrap;
    parameters["JrelNP0"] = m_ionicVariables.JrelNPcurrent;
    parameters["JrelCaMK0"] = m_ionicVariables.JrelCaMKcurrent;
    parameters["Nai0"] = m_ionicVariables.Nai;
    parameters["Nass0"] = m_ionicVariables.Nass;
    parameters["Ki0"] = m_ionicVariables.Ki;
    parameters["Kss0"] = m_ionicVariables.Kss;
    parameters["Cai0"] = m_ionicVariables.Cai;
    parameters["Cass0"] = m_ionicVariables.Cass;
    parameters["Cansr0"] = m_ionicVariables.Cansr;
    parameters["Cajsr0"] = m_ionicVariables.Cajsr;
    parameters["EAD"] = m_ionicVariables.ead;
    
    parameters["gKs"] = m_ionicVariables.gKs;
    parameters["gK1"] = m_ionicVariables.gK1;
    parameters["gNa_fast"] = m_ionicVariables.gNa_fast;
    parameters["gNa_late"] = m_ionicVariables.gNa_late;
    parameters["gto"] = m_ionicVariables.gto;
    parameters["gKr"] = m_ionicVariables.gKr;
    parameters["gCaL"] = m_ionicVariables.gCaL;
    parameters["gNaCa"] = m_ionicVariables.gNaCa;
    parameters["PCa"] = m_ionicVariables.PCa;// ADDED FOR SOBIE ARTICLE
    parameters["cellType"] = m_ionicVariables.cellType;// ADDED FOR SOBIE ARTICLE
    std::vector<double> sol_0;
    sol_0.resize(m_numODE);
    if(Vm_in_cvode) sol_0[m_numODE-1] = Vm0;
    
    sol_0[0] = parameters["m0"];
    sol_0[1] = parameters["hfast0"];
    sol_0[2] = parameters["hslow0"];
    sol_0[3] = parameters["j0"];
    sol_0[4] = parameters["hCaMKslow0"];
    sol_0[5] = parameters["jCaMK0"];
    sol_0[6] = parameters["mL0"];
    sol_0[7] = parameters["hL0"];
    sol_0[8] = parameters["hLCaMK0"];
    sol_0[9] = parameters["a0"];
    sol_0[10] = parameters["ifast0"];
    sol_0[11] = parameters["islow0"];
    sol_0[12] = parameters["aCaMK0"];
    sol_0[13] = parameters["iCaMKfast0"];
    sol_0[14] = parameters["iCaMKslow0"];
    sol_0[15] = parameters["d0"];
    sol_0[16] = parameters["ffast0"];
    sol_0[17] = parameters["fslow0"];
    sol_0[18] = parameters["fCafast0"];
    sol_0[19] = parameters["fCaslow0"];
    sol_0[20] = parameters["jCa0"];
    sol_0[21] = parameters["fCaMKfast0"];
    sol_0[22] = parameters["fCaCaMKfast0"];
    sol_0[23] = parameters["n0"];
    sol_0[24] = parameters["xrfast0"];
    sol_0[25] = parameters["xrslow0"];
    sol_0[26] = parameters["xs10"];
    sol_0[27] = parameters["xs20"];
    sol_0[28] = parameters["xk10"];
    sol_0[29] = parameters["CaMKtrap0"];
    sol_0[30] = parameters["JrelNP0"];
    sol_0[31] = parameters["JrelCaMK0"];
    sol_0[32] = parameters["Nai0"];
    sol_0[33] = parameters["Nass0"];
    sol_0[34] = parameters["Ki0"];
    sol_0[35] = parameters["Kss0"];
    sol_0[36] = parameters["Cai0"];
    sol_0[37] = parameters["Cass0"];
    sol_0[38] = parameters["Cansr0"];
    sol_0[39] = parameters["Cajsr0"];
    
    /* ADDED TO CHANGE INITIAL CONDITION */
    vector<double> heteroVec;
    if (cellHeterogeneity){
      readHeterogeneityField(heteroVec, pID);
    }
    ////////////////////////////////////////////
    if (cellHeterogeneity){
      /*
      for(int i=0 ; i<m_numVertexLocal ; i++){
        platformInt posglob = m_minOwnershipRange+i;
        std::vector<double> icCell = m_ionicVariables.m_InitialCond[getInitialCondition(heteroVec[posglob])];
        for(unsigned int iODE=0 ; iODE<m_vecSol.size() ; iODE++){
          VecSetValues(m_vecSol[iODE],1,&posglob,&icCell[iODE], INSERT_VALUES);
        }
      }
      */
      for(unsigned int iODE=0 ; iODE<m_vecSol.size() ; iODE++){
        VecSet(m_vecSol[iODE],sol_0[iODE]);
      }
    }
    else{
      for(unsigned int iODE=0 ; iODE<m_vecSol.size() ; iODE++){
        VecSet(m_vecSol[iODE],sol_0[iODE]);
      }
    }
    for(unsigned int iODE=0 ; iODE<m_vecSol.size() ; iODE++){
      VecAssemblyBegin(m_vecSol[iODE]);
      VecAssemblyEnd(m_vecSol[iODE]);
    }
    for(int i=0 ; i<m_numVertexLocal ; i++){
      platformInt posglob = m_minOwnershipRange+i;
      if (cellHeterogeneity){//if heterogeneous
        /*
        std::vector<double> icCell = m_ionicVariables.m_InitialCond[getInitialCondition(heteroVec[posglob])];
        for (size_t iODE=1 ; iODE<=m_vecSol.size() ; iODE++) {
          Ith(m_cvodeSol[i],iODE-1) = icCell[iODE-1];
          Ith(m_cvodeAbsTol[i],iODE-1) = m_atol;
        }
        */
        for (size_t iODE=1 ; iODE<=m_vecSol.size() ; iODE++) {        
          Ith(m_cvodeSol[i],iODE-1) = sol_0[iODE-1];
          Ith(m_cvodeAbsTol[i],iODE-1) = m_atol;
        }        
      }
      else {      
        for (size_t iODE=1 ; iODE<=m_vecSol.size() ; iODE++) {        
          Ith(m_cvodeSol[i],iODE-1) = sol_0[iODE-1];
          Ith(m_cvodeAbsTol[i],iODE-1) = m_atol;
        }
      }

      if (cellHeterogeneity){//if heterogeneous
        /*
        platformInt posglob = m_minOwnershipRange+i;
        std::vector<double> icCell = m_ionicVariables.m_InitialCond[getInitialCondition(heteroVec[posglob])];
        m_userData[i].m_Vm = icCell.back();
        */
        m_userData[i].m_Vm = Vm0;

        m_userData[i].m_cellFacto = getInitialCondition(heteroVec[posglob]);
      }
      else {
        m_userData[i].m_Vm = Vm0;
      }
      
      m_userData[i].m_VmInCVODE = m_VmInCVODE;
      m_userData[i].m_current_value.resize(m_listCurrent.size());
      m_userData[i].m_pos = m_minOwnershipRange+i;

      //For Strang splitting (VmInCVODE=true)
      m_userData[i].m_Cm = CmOHaraRudy;
      m_userData[i].m_dim0 = dim0;
      
      int flag = 0;
      flag = CVodeSetUserData(m_cvode_mem[i],&m_userData[i]);
      flag = CVodeInit(m_cvode_mem[i],IonicModel_ORd::M_f,0.0,m_cvodeSol[i]);
      flag = CVodeSVtolerances(m_cvode_mem[i],m_rtol,m_cvodeAbsTol[i]);
      
      if(m_VmInCVODE){
        flag = CVDense(m_cvode_mem[i],m_numODE);
      }
      else{
        flag = CVDense(m_cvode_mem[i],m_numODE-1);
      }
      flag = CVodeSetMaxStep(m_cvode_mem[i],m_maxTimeStep);
      flag = CVodeSetMaxOrd(m_cvode_mem[i],m_maxOrder);
      flag = CVodeSetMaxNumSteps(m_cvode_mem[i],m_maxNumStep);
    }
  }


  int IonicModel_ORd::M_f(realtype t, N_Vector y, N_Vector ydot, void *f_data) {

    /*
      Conductances:
      
      gKb = 0.003
      gKs = 0.0034
      gK1 = 0.1908
      gNa_fast = 75.0
      gNa_late = 0.0075
      gto = 0.02
      gKr = 0.046
      gNaCa = 0.0008
      gpCa = 0.0005
      
    */

    UserData* data = static_cast< UserData* >(f_data);
    int cellType = parameters["cellType"];

    double cellFacto = 1.0;
    
    if(data->m_cellHeterogeneity){
      cellType = data->m_cellType;
      const int idFacto = data->m_cellFacto;
      //cout << idFacto << endl;
      if (idFacto==1){
        cellFacto = 0.8;
      }
    }
    //cout << cellFacto << endl;
    double gKb = 0.003;
    double gKs = parameters["gKs"];
    double gK1 = parameters["gK1"];
    const double gNa_fast = parameters["gNa_fast"];
    double gNa_late = parameters["gNa_late"];
    double gto = parameters["gto"];
    double gKr = parameters["gKr"];
    double gNaCa = parameters["gNaCa"];
    const double gpCa = 0.0005;
    const double gCaL = parameters["gCaL"];

    double PCa = parameters["PCa"];//0.0001;//ADDED FOR SOBIE ARTICLE
    
    if (cellType==1){
      gNa_late *= 0.6;
      gto *= 4.0;
      gKr *= 1.3;
      gKs *= 1.4;
      gK1 *= 1.2;
      gNaCa *= 1.1;
      gKb *= 0.6;
      PCa *= 1.2;
      
    }
    else if (cellType==2){
      gto *= 4.0;
      gKr *= 0.8;
      gK1 *= 1.3;
      gNaCa *= 1.4;
      PCa *= 2.5;
      
    }

    gKr *= cellFacto;
    
    //UserData* data = static_cast< UserData* >(f_data);

    const double Na_o = 140.0;// mM
    const double Ca_o = 1.8;
    const double K_o = 5.4;
    
    double Vm;
    if(data->m_VmInCVODE){
      Vm = Ith(y,40);
    }
    else{
      Vm = data->m_Vm;
    }
    
    double iapp = 0.;
    //if(data->m_dim0){
    iapp = stimulation.value_Iapp(data->m_pos,t);
    //}
    
    double minf, taum;
    IonicModel_ORd::mgate(minf, taum, Vm);
    
    double hinf, tauhfast, tauhslow;
    IonicModel_ORd::hgate(hinf, tauhfast, tauhslow, Vm);
    
    double jinf, tauj;
    IonicModel_ORd::jgate(jinf, tauj, Vm);
    
    double hCaMKinf, tauhCaMKslow;
    IonicModel_ORd::hCaMKgate(hCaMKinf, tauhCaMKslow, Vm);
    
    double jCaMKinf, taujCaMKslow;
    IonicModel_ORd::jCaMKgate(jCaMKinf, taujCaMKslow, Vm);
    
    double mLinf, taumL;
    IonicModel_ORd::mLgate(mLinf, taumL, Vm);
    
    double hLinf, tauhL;
    IonicModel_ORd::hLgate(hLinf, tauhL, Vm);
    
    double hLCaMKinf, tauhLCaMK;
    IonicModel_ORd::hLCaMKgate(hLCaMKinf, tauhLCaMK, Vm);
    
    double ainf, taua;
    IonicModel_ORd::agate(ainf, taua, Vm);
    
    double iinf, tauifast, tauislow, aifast;
    IonicModel_ORd::igate(iinf, tauifast, tauislow, aifast, Vm);
    
    double aCaMKinf, tauaCaMK;
    IonicModel_ORd::aCaMKgate(aCaMKinf, tauaCaMK, Vm);
    
    double iCaMKinf, tauiCaMKfast, tauiCaMKslow, aiCaMKfast;
    IonicModel_ORd::iCaMKgate(iCaMKinf, tauiCaMKfast, tauiCaMKslow, aiCaMKfast, Vm);
    
    double dinf, taud;
    IonicModel_ORd::dgate(dinf, taud, Vm);
    
    double finf, tauffast, taufslow;
    IonicModel_ORd::fgate(finf, tauffast, taufslow, Vm);
    
    double fCainf, taufCafast, taufCaslow, afCafast;
    IonicModel_ORd::fCagate(fCainf, taufCafast, taufCaslow, afCafast, Vm);
    
    double jCainf, taujCa;
    IonicModel_ORd::jCagate(jCainf, taujCa, Vm);
    
    double fCaMKinf, taufCaMKfast;
    IonicModel_ORd::fCaMKgate(fCaMKinf, taufCaMKfast, Vm);
    
    double fCaCaMKinf, taufCaCaMKfast;
    IonicModel_ORd::fCaCaMKgate(fCaCaMKinf, taufCaCaMKfast, Vm);
    
    double xrinf, tauxrfast, tauxrslow, axrfast;
    IonicModel_ORd::xrgate(xrinf, tauxrfast, tauxrslow, axrfast, Vm);
    
    double xs1inf, tauxs1;
    IonicModel_ORd::xs1gate(xs1inf, tauxs1, Vm);
    
    double xs2inf, tauxs2;
    IonicModel_ORd::xs2gate(xs2inf, tauxs2, Vm);
    
    double xk1inf, tauxk1;
    IonicModel_ORd::xk1gate(xk1inf, tauxk1, Vm, K_o);
    
    const double L = 0.01;// cm
    const double rad = 0.0011;// cm
    const double vcell = 1000*3.14*rad*rad*L;// uL
    const double Ageo = 2*3.14*rad*rad+2*3.14*rad*L;// cm2
    const double Acap = 2*Ageo;// cm2
    const double vmyo = 0.68*vcell;// uL
    const double vnsr = 0.0552*vcell;
    const double vjsr = 0.0048*vcell;
    const double vss = 0.02*vcell;

    const double vf2rt = Vm*Frdy*Frdy/(RmOHaraRudy*Tcell);
    const double vfrt = Vm*Frdy/(RmOHaraRudy*Tcell);
    const double rtf = RmOHaraRudy*Tcell/Frdy;

    // Nernst Potentials
    const double ENa = rtf*log(Na_o/Ith(y, 32));
    const double EK = rtf*log(K_o/Ith(y, 34));
    const double EKs = rtf*log((K_o+0.01833*Na_o)/(Ith(y, 34)+0.01833*Ith(y, 32)));
    
    //------> CaMKactive///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const double alpha_CaMK = 0.05;// ms-1
    const double beta_CaMK = 0.00068;
    const double CaMK0 = 0.05;
    const double KmCaM = 0.0015;// mM
    
    const double facto = 1.0;
    const double CaMKbound = CaMK0*(1.-Ith(y, 29))/(1.+KmCaM/Ith(y, 37));// *facto pas mal pour Blebbistatin mais manque de variation sur le Cai!
    const double CaMKactive = facto*(CaMKbound + Ith(y, 29));
    Ith(ydot, 29) = alpha_CaMK*CaMKbound*(CaMKbound+Ith(y, 29))-beta_CaMK*Ith(y, 29);// CaMKtrap




    //INa fast current-----------------------------------------------------------------------------------------
    Ith(ydot, 0) = (minf-Ith(y, 0))/taum;// m_gate
    Ith(ydot, 1) = (hinf-Ith(y, 1))/tauhfast;// hfast_gate
    Ith(ydot, 2) = (hinf-Ith(y, 2))/tauhslow;// hslow_gate

    const double Ahf=0.99;
    const double Ahs=1.0-Ahf;
    const double h_gate = Ahf*Ith(y, 1)+Ahs*Ith(y, 2);// h_gate
    Ith(ydot, 3) = (jinf-Ith(y, 3))/tauj;// j_gate

    Ith(ydot, 4) = (hCaMKinf-Ith(y, 4))/tauhCaMKslow;// hCaMK_gate

    const double hCaMK_gate = Ahf*Ith(y, 1)+Ahs*Ith(y, 4);// hCaMK_gate

    Ith(ydot, 5) = (jCaMKinf-Ith(y,5))/taujCaMKslow;// jCaMK_gate


    const double KmCaMK = 0.15;
    const double O_INaCaMK = 1./(1.+(KmCaMK/CaMKactive));
    const double INa_fast = gNa_fast*(Vm-ENa)*pow(Ith(y, 0),3.)*((1.-O_INaCaMK)*h_gate*Ith(y, 3)+O_INaCaMK*hCaMK_gate*Ith(y,5));
    //const double INa_fast = 75.*(Vm-ENa)*pow(Ith(y, 0),3.)*((1.-O_INaCaMK)*h_gate*Ith(y, 3)+O_INaCaMK*hCaMK_gate*Ith(y,5));
   
    //---------------------------------------------------------------------------------------------------------

    //INa slow current-----------------------------------------------------------------------------------------
    Ith(ydot, 6) = (mLinf-Ith(y, 6))/taumL;// mL_gate
    Ith(ydot, 7) = (hLinf-Ith(y, 7))/tauhL;// hL_gate
    Ith(ydot, 8) = (hLCaMKinf-Ith(y, 8))/tauhLCaMK;// hLCaMK_gate
    
    const double INa_late = gNa_late*(Vm-ENa)*Ith(y, 6)*((1.-O_INaCaMK)*Ith(y, 7)+O_INaCaMK*Ith(y, 8));
    //const double INa_late = 0.0075*(Vm-ENa)*Ith(y, 6)*((1.-O_INaCaMK)*Ith(y, 7)+O_INaCaMK*Ith(y, 8));

    //---------------------------------------------------------------------------------------------------------

    const double INa = INa_fast + INa_late;

    //Ito current----------------------------------------------------------------------------------------------
    double delta_epi = 1.0;
    if (cellType==1){
      delta_epi =1.0-(0.95/(1.0+exp((Vm+70.0)/5.0)));
    }
    tauifast *= delta_epi;
    tauislow *= delta_epi;


    
    Ith(ydot, 9) = (ainf-Ith(y, 9))/taua;// a_gate
    Ith(ydot, 10) = (iinf-Ith(y, 10))/tauifast;// ifast_gate
    Ith(ydot, 11) = (iinf-Ith(y, 11))/tauislow;// islow_gate

    const double i_gate = aifast*Ith(y, 10)+(1.-aifast)*Ith(y, 11);// i_gate

    Ith(ydot, 12) = (aCaMKinf-Ith(y, 12))/tauaCaMK;// aCaMK_gate
    Ith(ydot, 13) = (iCaMKinf-Ith(y, 13))/tauiCaMKfast;// iCaMKfast_gate
    Ith(ydot, 14) = (iCaMKinf-Ith(y, 14))/tauiCaMKslow;// iCaMKslow_gate

    const double iCaMK_gate = aiCaMKfast*Ith(y, 13)+(1.-aiCaMKfast)*Ith(y, 14);// iCaMK_gate

    const double O_ItoCaMK = 1./(1.+(KmCaMK/CaMKactive));

    const double Ito = gto*(Vm-EK)*((1.-O_ItoCaMK)*Ith(y, 9)*i_gate+O_ItoCaMK*Ith(y, 12)*iCaMK_gate);
    //const double Ito = 0.02*(Vm-EK)*((1.-O_ItoCaMK)*Ith(y, 9)*i_gate+O_ItoCaMK*Ith(y, 12)*iCaMK_gate);

    //---------------------------------------------------------------------------------------------------------


    //ICaL current---------------------------------------------------------------------------------------------
    Ith(ydot, 15) = (dinf-Ith(y, 15))/taud;// d_gate
    Ith(ydot, 16) = (finf-Ith(y, 16))/tauffast;// ffast_gate
    Ith(ydot, 17) = (finf-Ith(y, 17))/taufslow;// fslow_gate

    const double Aff = 0.6;
    const double Afs = 1.-Aff;
    const double f_gate = Aff*Ith(y, 16)+Afs*Ith(y, 17);// f_gate

    Ith(ydot, 18) = (fCainf-Ith(y, 18))/taufCafast;// fCafast_gate
    Ith(ydot, 19) = (fCainf-Ith(y, 19))/taufCaslow;// fCaslow_gate


    const double afCaslow = 1.-afCafast;
    const double fCa_gate = afCafast*Ith(y, 18)+afCaslow*Ith(y, 19);// fCa_gate

    Ith(ydot, 20) = (jCainf-Ith(y, 20))/taujCa;// jCa_gate

    Ith(ydot, 21) = (fCaMKinf-Ith(y, 21))/taufCaMKfast;// fCaMKfast_gate
    const double fCaMKslow_gate = Ith(y, 17);// =fslow_gate

    const double fCaMK_gate = Aff*Ith(y, 21)+Afs*fCaMKslow_gate;// fCaMK_gate

    Ith(ydot, 22) = (fCaCaMKinf-Ith(y, 22))/taufCaCaMKfast;// fCaCaMKfast_gate
    const double fCaCaMKslow_gate = Ith(y, 19);// =fCaslow_gate

    const double fCaCaMK_gate = afCafast*Ith(y, 22)+(1.-afCafast)*fCaCaMKslow_gate;// fCaCaMK_gate
    
    const double Kmn = 0.002;
    const double kp2n = 1000.0;
    const double km2n = Ith(y, 20);

    const double alphan = 1./((kp2n/km2n)+pow(1.+(Kmn/Ith(y, 37)),4.));

    Ith(ydot, 23) = alphan*kp2n-Ith(y, 23)*km2n;// n_gate




    //[ NEW
    const double PhiCaL=4.0*vf2rt*(Ith(y, 37)*exp(2.0*vfrt)-0.341*Ca_o)/(exp(2.0*vfrt)-1.0);
    const double PhiCaNa=1.0*vf2rt*(0.75*Ith(y, 33)*exp(1.0*vfrt)-0.75*Na_o)/(exp(1.0*vfrt)-1.0);
    const double PhiCaK=1.0*vf2rt*(0.75*Ith(y, 35)*exp(1.0*vfrt)-0.75*K_o)/(exp(1.0*vfrt)-1.0);
    const double zca=2.0;
    //const double PCa=0.0001;
    const double PCap=1.1*PCa;
    const double PCaNa=0.00125*PCa;
    const double PCaK=3.574e-4*PCa;
    const double PCaNap=0.00125*PCap;
    const double PCaKp=3.574e-4*PCap;
    const double fICaLp=(1.0/(1.0+KmCaMK/CaMKactive));
    
    const double nCa = Ith(y, 23);
    const double jCa = Ith(y, 20);
    const double dCa = Ith(y, 15);
    const double fp = fCaMK_gate;
    const double fCap = fCaCaMK_gate;
 
    const double ICaL=gCaL*((1.0-fICaLp)*PCa*PhiCaL*dCa*(f_gate*(1.0-nCa)+jCa*fCa_gate*nCa)+fICaLp*PCap*PhiCaL*dCa*(fp*(1.0-nCa)+jCa*fCap*nCa));
    const double ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*dCa*(f_gate*(1.0-nCa)+jCa*fCa_gate*nCa)+fICaLp*PCaNap*PhiCaNa*dCa*(fp*(1.0-nCa)+jCa*fCap*nCa);
    const double ICaK=(1.0-fICaLp)*PCaK*PhiCaK*dCa*(f_gate*(1.0-nCa)+jCa*fCa_gate*nCa)+fICaLp*PCaKp*PhiCaK*dCa*(fp*(1.0-nCa)+jCa*fCap*nCa);
    
    //---------------------------------------------------------------------------------------------------------

    //IKr current----------------------------------------------------------------------------------------------
    Ith(ydot, 24) = (xrinf-Ith(y, 24))/tauxrfast;// xrfast_gate
    Ith(ydot, 25) = (xrinf-Ith(y, 25))/tauxrslow;// xrslow_gate

    const double axrslow = 1.-axrfast;
    const double xr_gate = axrfast*Ith(y, 24)+axrslow*Ith(y, 25);// xr_gate

    const double RKr = 1.0/(1.0+exp((Vm+55.0)/75.0))*1.0/(1.0+exp((Vm-10.0)/30.0));    
    const double IKr = parameters["EAD"]*gKr*sqrt(K_o/5.4)*xr_gate*RKr*(Vm-EK);////// NORMAL

    //---------------------------------------------------------------------------------------------------------


    //IKs current----------------------------------------------------------------------------------------------
    Ith(ydot, 26) = (xs1inf-Ith(y, 26))/tauxs1;// xs1_gate
    Ith(ydot, 27) = (xs2inf-Ith(y, 27))/tauxs2;// xs2_gate


    //const double gKs = 0.0034;
    const double IKs = gKs*(1.+0.6/(1.+pow(3.8e-5/Ith(y, 36),1.4)))*Ith(y, 26)*Ith(y, 27)*(Vm-EKs);
    //const double IKs = 0.0034*(1.+0.6/(1.+pow(3.8e-5/Ith(y, 36),1.4)))*Ith(y, 26)*Ith(y, 27)*(Vm-EKs);

    //---------------------------------------------------------------------------------------------------------

    //IK1 current----------------------------------------------------------------------------------------------
    Ith(ydot, 28) = (xk1inf-Ith(y, 28))/tauxk1;// xk1_gate

    const double RK1 = 1.0/(1.0+exp((Vm+105.8-2.6*K_o)/9.493));
    const double IK1 = gK1*sqrt(K_o)*Ith(y, 28)*RK1*(Vm-EK);
    //const double IK1 = 0.1908*sqrt(K_o)*Ith(y, 28)*RK1*(Vm-EK);

    //---------------------------------------------------------------------------------------------------------

    //INaCa current--------------------------------------------------------------------------------------------
    /*
    const double kNa1 = 15.0;// mM
    const double kNa2 = 5.0;
    const double kNa3 = 88.12;

    const double kasymm = 12.5;
    const double wNa = 6.0e4;// Hz
    const double wCa = 6.0e4;
    const double wNaCa = 5.0e3;

    const double kCaOn = 1.5e6;// mM/ms
    const double kCaOff = 5.0e3;// Hz
    const double qNa = 0.5224;
    const double qCa = 0.167;

    const double hCa = exp(qCa*vfrt);
    const double hNa = exp(qNa*vfrt);
    */
    const double KmCaAct = 150.0e-6;// mM-------ETAIT ICI AVANT
    
    const double kna1=15.0;
    const double kna2=5.0;
    const double kna3=88.12;
    const double kasymm=12.5;
    const double wna=6.0e4;
    const double wca=6.0e4;
    const double wnaca=5.0e3;
    const double kcaon=1.5e6;
    const double kcaoff=5.0e3;
    const double qna=0.5224;
    const double qca=0.1670;
    const double hca=exp(qca*vfrt);
    const double hna=exp(qna*vfrt);
    const double zna=1.0;
    
    double h1=1+Ith(y, 32)/kna3*(1+hna);
    double h2=(Ith(y, 32)*hna)/(kna3*h1);
    double h3=1.0/h1;
    double h4=1.0+Ith(y, 32)/kna1*(1+Ith(y, 32)/kna2);
    double h5=Ith(y, 32)*Ith(y, 32)/(h4*kna1*kna2);
    double h6=1.0/h4;
    double h7=1.0+Na_o/kna3*(1.0+1.0/hna);
    double h8=Na_o/(kna3*hna*h7);
    double h9=1.0/h7;
    double h10=kasymm+1.0+Na_o/kna1*(1.0+Na_o/kna2);
    double h11=Na_o*Na_o/(h10*kna1*kna2);
    double h12=1.0/h10;
    double k1=h12*Ca_o*kcaon;
    double k2=kcaoff;
    double k3p=h9*wca;
    double k3pp=h8*wnaca;
    double k3=k3p+k3pp;
    double k4p=h3*wca/hca;
    double k4pp=h2*wnaca;
    double k4=k4p+k4pp;
    double k5=kcaoff;
    double k6=h6*Ith(y, 36)*kcaon;
    double k7=h5*h2*wna;
    double k8=h8*h11*wna;
    double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
    double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
    double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
    double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
    double E1=x1/(x1+x2+x3+x4);
    double E2=x2/(x1+x2+x3+x4);
    double E3=x3/(x1+x2+x3+x4);
    double E4=x4/(x1+x2+x3+x4);
    double KmCaAct2 = 150.0e-6;
    double allo=1.0/(1.0+pow(KmCaAct2/Ith(y, 36),2.0));
    double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    double JncxCa=E2*k2-E1*k1;
    
    const double INaCa_i=0.8*gNaCa*allo*(zna*JncxNa+zca*JncxCa);
    
    h1=1+Ith(y, 33)/kna3*(1+hna);
    h2=(Ith(y, 33)*hna)/(kna3*h1);
    h3=1.0/h1;
    h4=1.0+Ith(y, 33)/kna1*(1+Ith(y, 33)/kna2);
    h5=Ith(y, 33)*Ith(y, 33)/(h4*kna1*kna2);
    h6=1.0/h4;
    h7=1.0+Na_o/kna3*(1.0+1.0/hna);
    h8=Na_o/(kna3*hna*h7);
    h9=1.0/h7;
    h10=kasymm+1.0+Na_o/kna1*(1+Na_o/kna2);
    h11=Na_o*Na_o/(h10*kna1*kna2);
    h12=1.0/h10;
    k1=h12*Ca_o*kcaon;
    k2=kcaoff;
    k3p=h9*wca;
    k3pp=h8*wnaca;
    k3=k3p+k3pp;
    k4p=h3*wca/hca;
    k4pp=h2*wnaca;
    k4=k4p+k4pp;
    k5=kcaoff;
    k6=h6*Ith(y, 37)*kcaon;
    k7=h5*h2*wna;
    k8=h8*h11*wna;
    x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
    x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
    x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
    x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
    E1=x1/(x1+x2+x3+x4);
    E2=x2/(x1+x2+x3+x4);
    E3=x3/(x1+x2+x3+x4);
    E4=x4/(x1+x2+x3+x4);
    KmCaAct2 = 150.0e-6;
    allo=1.0/(1.0+pow(KmCaAct2/Ith(y, 37),2.0));

    JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    JncxCa=E2*k2-E1*k1;
    
    const double INaCa_ss=0.2*gNaCa*allo*(zna*JncxNa+zca*JncxCa);
    
    const double INaCa=INaCa_i+INaCa_ss;
    
    //---------------------------------------------------------------------------------------------------------
    
    
    
    //INaK current---------------------------------------------------------------------------------------------
    
    const double k1p=949.5;
    const double k1m=182.4;
    const double k2p=687.2;
    const double k2m=39.4;
    k3p=1899.0;
    const double k3m=79300.0;
    k4p=639.0;
    const double k4m=40.0;
    const double Knai0=9.073;
    const double Knao0=27.78;
    const double delta=-0.1550;
    const double Knai=Knai0*exp(delta*vfrt/3.);
    const double Knao=Knao0*exp((1.-delta)*vfrt/3.);
    const double Kki=0.5;
    const double Kko=0.3582;
    const double MgADP=0.05;
    const double MgATP=9.8;
    const double Kmgatp=1.698e-7;
    const double H=1.0e-7;
    const double eP=4.2;
    const double Khp=1.698e-7;
    const double Knap=224.0;
    const double Kxkur=292.0;
    const double P=eP/(1.0+H/Khp+Ith(y, 32)/Knap+Ith(y, 34)/Kxkur);
    const double a1=(k1p*pow(Ith(y, 32)/Knai,3.0))/(pow(1.0+Ith(y, 32)/Knai,3.0)+pow(1.0+Ith(y, 34)/Kki,2.0)-1.0);
    const double b1=k1m*MgADP;
    const double a2=k2p;
    const double b2=(k2m*pow(Na_o/Knao,3.0))/(pow(1.0+Na_o/Knao,3.0)+pow(1.0+K_o/Kko,2.0)-1.0);
    const double a3=(k3p*pow(K_o/Kko,2.0))/(pow(1.0+Na_o/Knao,3.0)+pow(1.0+K_o/Kko,2.0)-1.0);
    const double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
    const double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
    const double b4=(k4m*pow(Ith(y, 34)/Kki,2.0))/(pow(1.0+Ith(y, 32)/Knai,3.0)+pow(1.0+Ith(y, 34)/Kki,2.0)-1.0);
    x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
    x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
    x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
    x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
    E1=x1/(x1+x2+x3+x4);
    E2=x2/(x1+x2+x3+x4);
    E3=x3/(x1+x2+x3+x4);
    E4=x4/(x1+x2+x3+x4);
    const double zk=1.0;
    const double JnakNa=3.0*(E1*a3-E2*b3);
    const double JnakK=2.0*(E4*b1-E3*a1);
    double Pnak = 30.;
    if (cellType==1){
      Pnak *= 0.9;
    }
    else if (cellType==2){
      Pnak *= 0.7;
    }
    
    const double INaK = Pnak*(zna*JnakNa+zk*JnakK);
    //---------------------------------------------------------------------------------------------------------


    //IpCa current---------------------------------------------------------------------------------------------

    //Formulation for INab, ICab, IKb and IpCa taken from Hund-Decker-Rudy model

    const double PNab = 3.75e-10;// cm/s
    //const double INab = PNab*zna*zna*vf2rt*(Ith(y, 32)*exp(zna*vfrt)-Na_o)/(exp(zna*vfrt)-1.0);
    const double INab = PNab*vf2rt*(Ith(y, 32)*exp(vfrt)-Na_o)/(exp(vfrt)-1.0);

    const double PCab = 2.5e-8;// cm/s
    const double ICab = PCab*zca*zca*vf2rt*(Ith(y, 36)*exp(2.0*vfrt)-0.341*Ca_o)/(exp(2.0*vfrt)-1.0);


    const double xKb = 1.0/(1.0+exp(-(Vm-14.48)/18.34));
    //const double IKb = gKb*xKb*(Vm-EK);
    
    const double IKb = gKb*xKb*(Vm-EK);
    
    //const double IpCa = gpCa*Ith(y, 36)/(0.0005+Ith(y, 36));
    const double IpCa = 0.0005*Ith(y, 36)/(0.0005+Ith(y, 36));

    //---------------------------------------------------------------------------------------------------------


    //Diffusion fluxes-----------------------------------------------------------------------------------------

    const double tau_diffNa = 2.0;// ms
    const double tau_diffK = 2.0;
    const double tau_difffCa = 0.2;

    const double JdiffNa = (Ith(y, 33)-Ith(y, 32))/tau_diffNa;
    const double JdiffCa = (Ith(y, 37)-Ith(y, 36))/tau_difffCa;
    const double JdiffK = (Ith(y, 35)-Ith(y, 34))/tau_diffK;

    //---------------------------------------------------------------------------------------------------------


    //SR Calcium Release Flux----------------------------------------------------------------------------------

    const double beta_t = 4.75;// ms
    const double arel = 0.5*beta_t;

    double JrelNPinf = -arel*ICaL/(1.+pow(1.5/Ith(y, 39),8.0));
    if (cellType==2){
      JrelNPinf *= 1.7;
    }
    
    double taurelNP = beta_t/(1.+0.0123/Ith(y, 39));
    /*
    if(taurelNP<0.005){
      taurelNP = 0.005;
    }
    */
    if(taurelNP<0.001){//SOBIE ARTICLE
      taurelNP = 0.001;
    }
    

    Ith(ydot, 30) = (JrelNPinf-Ith(y, 30))/taurelNP;// JrelNP_current

    const double beta_tCaMK = 1.25*beta_t;
    const double arelCaMK = 0.5*beta_tCaMK;
    double JrelCaMKinf = -ICaL*arelCaMK/(1.0+pow(1.5/Ith(y, 39),8.0));
    if (cellType==2){
      JrelCaMKinf *= 1.7;
    }
    
    double taurelCaMK = beta_tCaMK/(1.0+0.0123/Ith(y, 39));
    /*
    if(taurelCaMK<0.005){
      taurelCaMK = 0.005;
    }
    */
    if(taurelCaMK<0.001){//SOBIE ARTICLE
      taurelCaMK = 0.001;
    }
    
    Ith(ydot, 31) = (JrelCaMKinf-Ith(y, 31))/taurelCaMK;// JrelCaMK_current

    const double O_JrelCaMK = 1./(1.+(KmCaMK/CaMKactive));
    const double Jrel = (1.-O_JrelCaMK)*Ith(y, 30) + O_JrelCaMK*Ith(y, 31);//Bien pour la Blebbistatin mais bloque le recepteur a la Ryanodine et ne bloque pas de maniere directe la myosine

    //---------------------------------------------------------------------------------------------------------


    //Calcium Uptake (SERCA)-----------------------------------------------------------------------------------
    
    double Jupnp=0.004375*Ith(y, 36)/(Ith(y, 36)+0.00092);
    double Jupp=2.75*0.004375*Ith(y, 36)/(Ith(y, 36)+0.00092-0.00017);
    if (cellType==1){
      Jupnp *= 1.3;
      Jupp *= 1.3;
    }

    
    const double fJupp=(1.0/(1.0+KmCaMK/CaMKactive));
    const double Jleak=0.0039375*Ith(y, 38)/15.0;
    const double Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;
    
    //---------------------------------------------------------------------------------------------------------

    //Calcium Translocation (NSR->JSR)-------------------------------------------------------------------------
    const double tautr = 100.0;// ms
    const double Jtr = (Ith(y, 38)-Ith(y, 39))/tautr;
    //---------------------------------------------------------------------------------------------------------


    //Concentrations-------------------------------------------------------------------------------------------

    double CMDN = 0.05;// mM
    if (cellType==1){
      CMDN *= 1.3;
    }
    
    const double KmCMDN = 0.00238;// mM

    const double TRPN = 0.07;
    const double KmTRPN = 0.0005;

    const double BSR = 0.047;
    const double KmBSR = 0.00087;

    const double BSL = 1.124;
    const double KmBSL = 0.0087;

    const double CSQN = 10.0;
    const double KmCSQN = 0.8;



    Ith(ydot, 32) = -(INa_fast+INa_late+3.0*INaCa_i+3.0*INaK+INab)*Acap/(Frdy*vmyo)+JdiffNa*vss/vmyo;// Na_i
    Ith(ydot, 33) = -(ICaNa+3.0*INaCa_ss)*Acap/(Frdy*vss)-JdiffNa;// Na_ss

    Ith(ydot, 34) = -(Ito+IKr+IKs+IK1+IKb+iapp-2.*INaK)*(Acap/(Frdy*vmyo))+JdiffK*vss/vmyo;// K_i
    Ith(ydot, 35) = -ICaK*(Acap/(Frdy*vss))-JdiffK;// K_ss
    
    const double beta_Cai = 1. / (1.+(CMDN*KmCMDN/pow(KmCMDN+Ith(y, 36),2.))+(TRPN*KmTRPN/pow(KmTRPN+Ith(y, 36),2.)));
    Ith(ydot, 36) = beta_Cai*(-(IpCa+ICab-2.*INaCa_i)*Acap/(2.0*Frdy*vmyo)-Jup*vnsr/vmyo+JdiffCa*vss/vmyo);//Ca_i


    const double beta_Cass = 1. / (1.+(BSR*KmBSR/pow(KmBSR+Ith(y, 37),2.))+(BSL*KmBSL/pow(KmBSL+Ith(y, 37),2.)));
    Ith(ydot, 37) = beta_Cass*(-(ICaL-2.*INaCa_ss)*Acap/(2.0*Frdy*vss)+Jrel*vjsr/vss-JdiffCa);// Ca_ss

    Ith(ydot, 38) = Jup-Jtr*vjsr/vnsr;// Ca_nsr

    const double beta_Cajsr = 1./(1.+(CSQN*KmCSQN/pow(KmCSQN+Ith(y, 39),2.)));
    Ith(ydot, 39) = beta_Cajsr*(Jtr-Jrel);// Ca_jsr

    //---------------------------------------------------------------------------------------------------------

    double value_I_ion = -(INa_fast+INa_late+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+iapp);// /CmOHaraRudy

   
    if(data->m_VmInCVODE) Ith(ydot, 40) = value_I_ion;
    data->m_value_I_ion = value_I_ion;
    
    
    data->m_current_value[0] = INa;
    data->m_current_value[1] = Ito;
    data->m_current_value[2] = ICaL;
    data->m_current_value[3] = ICaNa;
    data->m_current_value[4] = ICaK;
    data->m_current_value[5] = IKr;
    data->m_current_value[6] = IKs;
    data->m_current_value[7] = IK1;
    data->m_current_value[8] = INaCa_i;
    data->m_current_value[9] = INaCa_ss;
    data->m_current_value[10] = INaCa;
    data->m_current_value[11] = INaK;
    data->m_current_value[12] = INab;
    data->m_current_value[13] = ICab;
    data->m_current_value[14] = IKb;
    data->m_current_value[15] = IpCa;
    data->m_current_value[16] = JdiffNa;
    data->m_current_value[17] = JdiffK;
    data->m_current_value[18] = JdiffCa;
    data->m_current_value[19] = Ith(y, 30);
    data->m_current_value[20] = Ith(y, 31);
    data->m_current_value[21] = Jrel;
    data->m_current_value[22] = Jupnp;
    data->m_current_value[23] = Jupp;
    data->m_current_value[24] = Jup;
    data->m_current_value[25] = Jtr;
    
    return 0;
    
  }

#pragma endregion
  
  //__________________________________________________________________________________________________________________________________________________________________________
  //__________________________________________________________________________________________________________________________________________________________________________

  

  

  
  //__________________________________________________________________________________________________________________________________________________________________________
  //__________________________________________________________________________________________________________________________________________________________________________


  //*******//
  // GATES //
  //*******//

  
#pragma region ORd gates
  
  void IonicModel_ORd::mgate(double& minf, double& taum, double u) {
    minf = 1.0/(1.0+exp((-(u+39.57))/9.871));
    taum = 1.0/(6.765*exp((u+11.64)/34.77)+8.552*exp(-(u+77.42)/5.955));
  }
  
  void IonicModel_ORd::hgate(double& hinf, double& tauhfast, double& tauhslow, double u) {
    hinf = 1.0/(1+exp((u+82.90)/6.086));
    tauhfast = 1.0/(1.432e-5*exp(-(u+1.196)/6.285)+6.149*exp((u+0.5096)/20.27));
    tauhslow = 1.0/(0.009794*exp(-(u+17.95)/28.05)+0.3343*exp((u+5.730)/56.66));
  }
  
  void IonicModel_ORd::jgate(double& jinf, double& tauj, double u) {
    jinf = 1.0/(1+exp((u+82.90)/6.086));
    tauj = 2.038+1.0/(0.02136*exp(-(u+100.6)/8.281)+0.3052*exp((u+0.9941)/38.45));
  }
  
  void IonicModel_ORd::hCaMKgate(double& hCaMKinf, double& tauhCaMKslow, double u) {
    hCaMKinf = 1.0/(1+exp((u+89.1)/6.086));
    tauhCaMKslow = 3.0*(1.0/(0.009794*exp(-(u+17.95)/28.05)+0.3343*exp((u+5.730)/56.66)));
  }
  
  void IonicModel_ORd::jCaMKgate(double& jCaMKinf, double& taujCaMK, double u) {
    jCaMKinf = 1.0/(1+exp((u+82.90)/6.086));
    taujCaMK = 1.46*(2.038+1.0/(0.02136*exp(-(u+100.6)/8.281)+0.3052*exp((u+0.9941)/38.45)));
  }
  
  void IonicModel_ORd::mLgate(double& mLinf, double& taumL, double u) {
    mLinf = 1.0/(1.0+exp((-(u+42.85))/5.264));
    taumL = 1.0/(6.765*exp((u+11.64)/34.77)+8.552*exp(-(u+77.42)/5.955));
  }
  
  void IonicModel_ORd::hLgate(double& hLinf, double& tauhL, double u) {
    hLinf = 1.0/(1.0+exp((u+87.61)/7.488));
    tauhL = 200.0;
  }
  
  void IonicModel_ORd::hLCaMKgate(double& hLCaMKinf, double& tauhLCaMK, double u) {
    hLCaMKinf = 1.0/(1.0+exp((u+93.81)/7.488));
    tauhLCaMK = 3.0*200.0;
  }
  
  void IonicModel_ORd::agate(double& ainf, double& taua, double u) {
    ainf = 1.0/(1.0+exp((-(u-14.34))/14.82));
    taua = 1.0515/(1.0/(1.2089*(1.0+exp(-(u-18.4099)/29.3814)))+3.5/(1.0+exp((u+100.0)/29.3814)));
  }
  
  void IonicModel_ORd::igate(double& iinf, double& tauifast, double& tauislow, double& aifast, double u) {
    iinf = 1.0/(1.0+exp((u+43.94)/5.711));
    tauifast = 4.562+1/(0.3933*exp((-(u+100.0))/100.0)+0.08004*exp((u+50.0)/16.59));
    tauislow = 23.62+1/(0.001416*exp((-(u+96.52))/59.05)+1.780e-8*exp((u+114.1)/8.079));
    aifast = 1.0/(1.0+exp((u-213.6)/151.2));
  }
  
  void IonicModel_ORd::aCaMKgate(double& aCaMKinf, double& tauaCaMK, double u) {
    aCaMKinf = 1.0/(1.0+exp((-(u-24.34))/14.82));
    tauaCaMK = 1.0515/(1.0/(1.2089*(1.0+exp(-(u-18.4099)/29.3814)))+3.5/(1.0+exp((u+100.0)/29.3814)));
  }
  
  void IonicModel_ORd::iCaMKgate(double& iCaMKinf, double& tauiCaMKfast, double& tauiCaMKslow, double& aiCaMKfast, double u) {
    iCaMKinf = 1./(1.+exp((u+43.94)/5.711));
    
    double dummy_tauiCaMKfast = 4.562+1/(0.3933*exp((-(u+100.0))/100.0)+0.08004*exp((u+50.0)/16.59));
    double dummy_tauiCaMKslow = 23.62+1/(0.001416*exp((-(u+96.52))/59.05)+1.780e-8*exp((u+114.1)/8.079));
    
    double deltaDevelop = 1.354+1.0e-4/(exp((u-167.4)/15.89)+exp(-(u-12.23)/0.2154));
    double deltaRecover = 1.0-0.5/(1.0+exp((u+70.0)/20.0));
    
    tauiCaMKfast = dummy_tauiCaMKfast*deltaDevelop*deltaRecover;
    tauiCaMKslow = dummy_tauiCaMKslow*deltaDevelop*deltaRecover;
    
    aiCaMKfast = 1.0/(1.0+exp((u-213.6)/151.2));
    
  }
  
  void IonicModel_ORd::dgate(double& dinf, double& taud, double u) {
    dinf = 1.0/(1.0+exp((-(u+3.940))/4.230));
    taud = 0.6+1.0/(exp(-0.05*(u+6.0))+exp(0.09*(u+14.0)));
  }
  
  void IonicModel_ORd::fgate(double& finf, double& tauffast, double& taufslow, double u) {
    finf = 1.0/(1.0+exp((u+19.58)/3.696));
    tauffast = 7.0+1.0/(0.0045*exp(-(u+20.0)/10.0)+0.0045*exp((u+20.0)/10.0));
    taufslow = 1000.0+1.0/(0.000035*exp(-(u+5.0)/4.0)+0.000035*exp((u+5.0)/6.0));
  }
  
  void IonicModel_ORd::fCagate(double& fCainf, double& taufCafast, double& taufCaslow, double& afCafast, double u) {
    fCainf = 1.0/(1.0+exp((u+19.58)/3.696));
    taufCafast = 7.0+1.0/(0.04*exp(-(u-4.0)/7.0)+0.04*exp((u-4.0)/7.0));
    taufCaslow = 100.0+1.0/(0.00012*exp(-u/3.0)+0.00012*exp(u/7.0));
    afCafast = 0.3+0.6/(1.0+exp((u-10.0)/10.0));
  }
  
  void IonicModel_ORd::jCagate(double& jCainf, double& taujCa, double u) {
    jCainf = 1.0/(1.0+exp((u+19.58)/3.696));
    taujCa = 75.0;
  }
  
  void IonicModel_ORd::fCaMKgate(double& fCaMKinf, double& taufCaMKfast, double u) {
    fCaMKinf = 1.0/(1.0+exp((u+19.58)/3.696));
    taufCaMKfast = 2.5*(7.0+1.0/(0.0045*exp(-(u+20.0)/10.0)+0.0045*exp((u+20.0)/10.0)));
  }
  
  void IonicModel_ORd::fCaCaMKgate(double& fCaCaMKinf, double& taufCaCaMKfast, double u) {
    fCaCaMKinf = 1.0/(1.0+exp((u+19.58)/3.696));
    taufCaCaMKfast = 2.5*(7.0+1.0/(0.04*exp(-(u-4.0)/7.0)+0.04*exp((u-4.0)/7.0)));
  }
  
  void IonicModel_ORd::xrgate(double& xrinf, double& tauxrfast, double& tauxrslow, double& axrfast, double u) {
    xrinf = 1.0/(1.0+exp((-(u+8.337))/6.789));
    tauxrfast = 12.98+1.0/(0.3652*exp((u-31.66)/3.869)+4.123e-5*exp((-(u-47.78))/20.38));
    tauxrslow = 1.865+1.0/(0.06629*exp((u-34.70)/7.355)+1.128e-5*exp((-(u-29.74))/25.94));
    axrfast = 1.0/(1.0+exp((u+54.81)/38.21));
  }
  
  void IonicModel_ORd::xs1gate(double& xs1inf, double& tauxs1, double u) {
    xs1inf = 1.0/(1.0+exp((-(u+11.60))/8.932));
    tauxs1 = 817.3+1.0/(2.326e-4*exp((u+48.28)/17.80)+0.001292*exp((-(u+210.0))/230.0));
  }
  
  void IonicModel_ORd::xs2gate(double& xs2inf, double& tauxs2, double u) {
    xs2inf = 1.0/(1.0+exp((-(u+11.60))/8.932));
    tauxs2 = 1.0/(0.01*exp((u-50.0)/20.0)+0.0193*exp((-(u+66.54))/31.0));
  }
  
  void IonicModel_ORd::xk1gate(double& xk1inf, double& tauxk1, double u, double K) {
    xk1inf = 1.0/(1.0+exp(-(u+2.5538*K+144.59)/(1.5692*K+3.8115)));
    tauxk1 = 122.2/(exp((-(u+127.2))/20.36)+exp((u+236.8)/69.33));
  }
  
#pragma endregion
  

  //*******//
  //*******//
  
  
  //__________________________________________________________________________________________________________________________________________________________________________
  //__________________________________________________________________________________________________________________________________________________________________________
  
  //__________________________________________________________________________________________________________________________________________________________________________
  //__________________________________________________________________________________________________________________________________________________________________________
  
  }
