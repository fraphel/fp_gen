#include "stdafx.h"


#include "bidomainproblem.hpp"

namespace cardioxcomp {


  //#############//
  // CONSTRUCTOR //
  //#############//

#pragma region VOID CONSTRUCTOR
  Bidomainproblem::Bidomainproblem(){
  }
#pragma endregion
  
#pragma region LUA VARIABLES CONSTRUCTOR
  Bidomainproblem::Bidomainproblem(std::string meshfile,int numProc,int rankProc, const LuaFile& parameters, const LuaVariables& luavar, const string procId):
    Problem(meshfile,2,2,3,numProc,rankProc, parameters.isHeterogeneity, procId, parameters.cellHeterogeneity){
    
    m_parameters = parameters;
    m_variables = luavar;
    
    //Output
    m_writeODE = m_parameters.writeEDOs;
    m_writeCurrent = m_parameters.writeCurrents;
    m_procId = procId;

    m_physicalvariables.clear();
    m_physicalvariables.push_back("TransmembranarPotential");
    m_physicalvariables.push_back("FieldPotential");
    
    //FrequencyWriteSolution
    m_freq_writeSolution = m_parameters.freq_writeSolution;
    
    //Parameters for CVODE (relative/absolute tolerances, options:order, number of steps)
    m_rtol = m_parameters.rtol;
    m_atol = m_parameters.atol;
    m_maxOrder = m_parameters.maxOrder;
    m_maxNumStep = m_parameters.maxNumStep;
    m_maxStep = m_parameters.maxStep;
    
    //Parameters for the Solver (Ax=b)
    m_edpReltol = m_parameters.relativeTolerance;
    m_edpAbstol = m_parameters.absoluteTolerance;
    m_edpMaxiter = m_parameters.maxIteration;
    m_solver = m_parameters.solver;
    m_preconditioner = m_parameters.preconditioner;
    m_preconditionerOption = m_parameters.preconditionerOption;
    
    //Stimulation
    if(m_parameters.manual == true){
      m_stim = Stimulation(p_mesh,m_parameters.manual,2,m_parameters.starts,m_parameters.ends,m_parameters.iapps,m_parameters.focus);
    }
    else{
      m_stim = Stimulation(p_mesh,m_parameters.manual,2,m_parameters.start,m_parameters.iapp,m_parameters.duration,m_parameters.repetition,m_parameters.focus);
    }
    
    //Heterogeneity
    m_isHeterogeneity = m_parameters.isHeterogeneity;
    if(m_isHeterogeneity){
      m_heterogeneity = Heterogeneity(p_mesh,m_parameters.posHeterogeneity);

	  m_isAtrialVentricular = m_parameters.isAtrialVentricular;
	  m_isRescaling = m_parameters.isRescaling;
    }

    //Initial Vm heterogeneity
    if(m_ionicmodel == "OHaraRudy"){
      m_heteroInitial.clear();
      if (m_parameters.isHeterogeneity){
        vector<double> heteroVec;
        readHeterogeneityField(heteroVec, procId);
        for (int i=0 ; i<m_localsize/m_numComp ; i++) {
          platformInt pos = i;//i+m_minOwnership/m_numComp;
          platformInt idof = pos*p_dof->numComp();
        AOApplicationToPetsc(p_dof->ao(), 1, &idof);
        std::vector<double> icCell = luavar.m_InitialCond[getInitialCondition(heteroVec[pos])];
        VecSetValue(p_vec[0], i, icCell.back(), INSERT_VALUES); // initial values in p_vec[0]
        }
      }
    }
    
    //CellHeterogeneity
    m_cellHeterogeneity = m_parameters.cellHeterogeneity;
    
    //Parameters for the ionic model
    m_ionicmodel = m_parameters.ionicmodel;
    if(m_ionicmodel == "Paci"){
      if(m_parameters.cells == "Atrial"){
        m_atrial = true;
      }
      else if(m_parameters.cells == "Ventricular"){
        m_atrial = false;
      }
    else{
      std::cout << "WARNING: " << m_parameters.cells << " is a wrong parameter....: Atrial or Ventricular" << std::endl;
    }
    
    m_spontaneous = m_parameters.spontaneous;
    }
    
    //Conductivities
    m_sigmai = m_parameters.sigma_i;//intracellular
    m_sigmae = m_parameters.sigma_e;//extracellular
    
    m_Am = m_parameters.Am;
    m_Cm = m_parameters.Cm;
    
    //Coefficients for the PDE (problem.cpp)
    m_coefftime = m_Am*m_Cm;
    m_coeffsource = m_Am;
    
    m_dt = m_parameters.dt;
    m_Nmax = m_parameters.nmax;
    m_Tmax = m_parameters.tmax;
    m_orderbdf = m_parameters.VmOrderBdf;
    m_name = "bidomain";
    
    if(m_parameters.splitting){
      m_theta = parameters.theta_split;
    }
    
    //Initials values for Vm & Ue
    m_physicalInitialValues.clear();
    //m_physicalInitialValues.push_back(m_parameters.Vm_init);
    //m_physicalInitialValues.push_back(m_parameters.Ue_init);
   
    if(!m_isHeterogeneity){
      m_physicalInitialValues.push_back(m_parameters.Vm_init);
      m_physicalInitialValues.push_back(m_parameters.Ue_init);
    }
    else{
      if(!m_isAtrialVentricular){
        m_physicalInitialValues.push_back(m_parameters.Vm_init);
        m_physicalInitialValues.push_back(m_parameters.Ue_init);
      }
      else{
        m_physicalInitialValues.push_back(m_parameters.Vm_init);
        m_physicalInitialValues.push_back(m_parameters.Ue_init);
      }
    }


    m_nbEl = m_parameters.nbEl;
    SolutionElectrodes.resize(m_nbEl);
    SolutionVm.resize(m_nbEl);
    
    m_Uemean.resize(m_nbEl);
    m_Vmmean.resize(m_nbEl);
    m_UmesVec.resize(m_nbEl);
    m_VmVec.resize(m_nbEl);
    
    for(int ie=0 ; ie<m_nbEl ; ie++){
      m_Uemean[ie] = 0.;
      m_Vmmean[ie] = m_physicalInitialValues[0];
    }
    //m_CaimeanWell = 1e-4;
    
    if(m_ionicmodel=="FHN" || m_ionicmodel=="MV" || m_ionicmodel=="PCBE"){
      m_vmin = m_parameters.Vmin;
      m_vmax = m_parameters.Vmax;
      m_dv = m_vmax-m_vmin;
    }
    
    
    
    // Initialize electrodes parameters
    //_________________________________________________
    
    initializeElectrodeParameters(m_nbEl, m_parameters.Ri, m_parameters.Rel, m_parameters.Cel, m_parameters.diametres);
    m_electrodes.initializeElectrodes(m_nbEl,m_dt,m_Cel,m_Ri,m_Rel,m_diametres);
    //writeElectrodesSolution_headers();
    //writeElectrodesSolution();
    
    //_________________________________________________
    //_________________________________________________
    
    
    //For boundary condition
    m_labelsBC = m_parameters.labels;
    m_variablesBC = m_parameters.variables;
    m_valuesBC = m_parameters.values;
    
    //Solver of ionic model for Paci
    if(m_ionicmodel == "Paci"){
      m_isCVOdeSol = m_parameters.CVodeSolver;
    }
    
  }
#pragma endregion
  
  //#############//
  //#############//


  //#########################//
  // COMPUTE MATRICES/VECTOR //
  //#########################//
  
#pragma region COMPUTE MAT-VEC
  
  /*
    WARNING: in Ax=b, x=(Vm0,Ue0,Vm1,Ue1,....,Vmn,Uen)
    In local (for ionic models): (Vm0,Vm1,....Vmn)
  */

  void Bidomainproblem::computeMatVec(){
    int iglob,jglob;
    Mesh& localMesh = *p_localMesh;
    Mat& matrix = p_mat[0];  // p_mat is a table of Mat
    Mat& matrixMass = p_mat[1];
    Vec& sol = p_vec[0];  // p_vec is a table of Vec
    Vec& rhs = p_vec[1];
    Vec& fct = p_vec[2];
    
    double* matElemBidomain = p_matElem[0].val;  // p_matElem is of type MatElem (class defined in matElem.hpp)
    double* matElemMass = p_matElem[1].val;
    double* vecElem = p_vecElem[0].val;
    PetscInt* m_rowElem = new PetscInt[3*m_numComp];
    PetscInt* m_colElem = new PetscInt[3*m_numComp];

    MatElem mat;
    mat.init(m_numComp);

    for(int it=0;it<localMesh.nt;it++){  // cycle on the triangles
      const Triangle& tria = localMesh[it];
      for(int i=0;i<3;i++){  // cycle on the nodes of the triangle
        Point Gi = tria.H(i);

        for(int icomp=0;icomp<m_numComp;icomp++){  // cycle on the nodes globally
          iglob=localMesh(it,i) * m_numComp + icomp;
          m_rowElem[i*m_numComp+icomp] = iglob;
          for(int j=0;j<3;j++){  // cycle on the nodes of the triangle
            Point Gj = tria.H(j);

            for(int jcomp=0;jcomp<m_numComp;jcomp++){
              jglob=localMesh(it,j) * m_numComp + jcomp;
              m_colElem[j*m_numComp+jcomp] = jglob;

              int index = mat.index(icomp,jcomp,i,j);
              const double K_tria = (Gi,Gj)*tria.area;  // ?? should be the matrix for the gradient

              
              if(icomp == 0 && jcomp == 0){ // first diagonal block of the matrix in the coupled method
                const double coef = m_bdfEdp.coeffDeriv0()/m_dt;  // derivative coefficient from bdf method: with BE (bdf=1) is 1/m_dt
                  
                if(m_parameters.splitting){// if splitting method
                  matElemBidomain[index] = (m_Am*m_Cm*coef*(1. + (i==j)) /12.)*tria.area + m_theta*m_sigmai*K_tria;//Am*Cm*coeff*V(t+dt)+sigmai*grad(Vm)
                  matElemMass[index] = (1. + (i==j)) /12. * tria.area - m_dt*(1.-m_theta)*m_sigmai*K_tria/(m_Am*m_Cm);//For rhs: Vn-1 (Am*Cm to be consistant with the RHS (problem.cpp) )
                }
                else{
                  matElemBidomain[index] = (m_Am*m_Cm*coef*(1. + (i==j)) /12.)*tria.area + m_sigmai*K_tria;//Am*Cm*coeff*V(t+dt)+sigmai*grad(Vm)
                  matElemMass[index] = (1. + (i==j)) /12. * tria.area;//For rhs: Vn-1
                }
                
              }
              else if(icomp == 1 && jcomp == 1){ // second diagonal block of the matrix in the coupled method
                
                if(m_parameters.splitting){
                  matElemBidomain[index] = (m_sigmai + m_sigmae)*K_tria/m_theta;//(sigmai+sigmae)*grad(Ue)/theta
                }
                else{
                  matElemBidomain[index] = (m_sigmai + m_sigmae)*K_tria;//(sigmai+sigmae)*grad(Ue)
                }
                
              }
              else{ // blocks out of the diagonal
                matElemBidomain[index] = m_sigmai*K_tria;//sigmai*grad(Ue) and sigmai*grad(Vm) ---------------------------> B et transp(B)
                
                if(m_parameters.splitting==true && icomp==1){//Splitting
                  matElemMass[index] = -(1.-m_theta)*m_dt*m_sigmai*K_tria/(m_Am*m_Cm*m_theta);//For rhs: Vn-1 (Am*Cm to be consistant with the RHS (problem.cpp) )
                }
                
              }
            }
          }
          vecElem[i*m_numComp + icomp]= 0.;//source term: Iion, Iapp and Iel (calculated in Source function)
        }
      }
      AOApplicationToPetsc(p_dof->ao(),3*m_numComp,m_rowElem);
      AOApplicationToPetsc(p_dof->ao(),3*m_numComp,m_colElem);
      MatSetValues(matrix,3*m_numComp,m_rowElem,3*m_numComp,m_colElem,matElemBidomain,ADD_VALUES);
      MatSetValues(matrixMass,3*m_numComp,m_rowElem,3*m_numComp,m_colElem,matElemMass,ADD_VALUES);
      VecSetValues(fct,3*m_numComp,m_rowElem,vecElem, INSERT_VALUES);
    }
  }
  
#pragma endregion

  //#########################//
  //#########################//


  //###########//
  // DIRICHLET //
  //###########//

#pragma region DIRICHLET BC
	/*
	Applies Dirichlet boundary condition defined in the config.lua
	*/
  
  void Bidomainproblem::dirichletBC(){

    Mat& matrix = p_mat[0];
    Vec& rhs = p_vec[1];

    for(unsigned int ilab=0 ; ilab<m_labelsBC.size() ; ilab++){
      int var;
      if(m_variablesBC[ilab] == "extracellular"){
        var = 1;
      }
      else if(m_variablesBC[ilab] == "transmembranar"){
        var = 0;
      }
      else{
        std::cout << "ERROR in luafile: " << m_variablesBC[ilab] << " is a wrong variable name.....: extracellular or transmembranar" << std::endl;
        exit(0);
      }
      applyDirichletBC(matrix,rhs,m_valuesBC[ilab],m_labelsBC[ilab],*p_mesh,*p_dof,m_rankProc,var);
    }

  }

#pragma endregion

  //###########//
  //###########//


  //############//
  // DESTRUCTOR //
  //############//

#pragma region DESTRUCTOR
  Bidomainproblem::~Bidomainproblem(){

    delete m_ionicModel;
    
    for(unsigned int ie=0 ; ie<m_list_ODE.size() ; ie++){
      VecDestroy(&m_vecODE[ie]);
      VecDestroy(&m_vecSeqODE[ie]);
    }
    
    for(unsigned int ic=0 ; ic<m_list_Current.size() ; ic++){
      VecDestroy(&m_vecCurrent[ic]);
      VecDestroy(&m_vecSeqCurrent[ic]);
    }
    
    VecDestroy(&m_ion);
    
  }
#pragma endregion

  //############//
  //############//


  //########################//
  // Initialize ionic model //
  //########################//

#pragma region INITIALIZE PROBLEM

  void Bidomainproblem::initializeProblem(const Vec& uExtrap){
    
    //Create ionic model
    //______________________________
    
    if(m_ionicmodel == "OHaraRudy"){
      m_ionicModel = new IonicModel_ORd(0.,m_dt,p_mesh->nv,m_localsize/m_numComp,m_rtol,m_atol,m_maxOrder,m_maxNumStep,m_maxStep,m_parameters.splitting,m_physicalInitialValues[0], false, m_variables, m_procId,m_cellHeterogeneity);
    }
    else{
      std::cout << "Error config.lua.....IonicModel.Model = OHaraRudy only for this version" << std::endl;
      exit(0);
    }
    
    //setting the heterogeneity
    //m_ionicModel->setIsHeterogeneity(m_isHeterogeneity);

    if(m_ionicmodel != "Paci"){
      m_ionicModel->setIsHeterogeneity(m_isHeterogeneity);
    }


    if(m_isHeterogeneity){

      m_ionicModel->setHeterogeneity(m_heterogeneity);

      if(m_ionicmodel == "Paci"){
        m_ionicModel->setIsAtrialVentricular(m_isAtrialVentricular);
        m_ionicModel->setIsRescaling(m_isRescaling);
      }
    }
    
    //______________________________
    //______________________________

    //From global to local numerotation

    Vec Vinit;
    VecCreate(PETSC_COMM_WORLD,&Vinit);
    VecSetSizes(Vinit,m_localsize/m_numComp,p_mesh->nv);  // Vinit has as size the number of nodes
    VecSetFromOptions(Vinit);
    
    for(int i=0 ; i<m_localsize/m_numComp ; i++){

      platformInt idofloc = i*m_numComp;
      platformInt idof;
#ifdef PETSC_3_6
      ISLocalToGlobalMappingApply(p_dof->locToGlob(),1,&idofloc,&idof);
#else
      ISLocalToGlobalMappingApply(p_dof->return_mappingNodes(),1,&idofloc,&idof);
#endif
      AOApplicationToPetsc(p_dof->ao(), 1, &idof);//ordering application
      
      double value;
      VecGetValues(uExtrap,1,&idof,&value);   // extrapolated potential inside the variable value

      platformInt pos = i+m_minOwnership/m_numComp;
      
      // change the initial Vm for the heterogeneity part
      if(m_isHeterogeneity){
        if(m_heterogeneity.isInHetero(pos)){
          if(m_isAtrialVentricular){
            value = -0.068733823452164;
          }
          if(m_isRescaling){
            value = value*0.8 + m_parameters.Vm_init*0.2;
          }
        }
      }
      
      VecSetValues(Vinit,1,&pos,&value, INSERT_VALUES); // potential inside Vinit

    }
    
    VecAssemblyBegin(Vinit);
    VecAssemblyEnd(Vinit);
    
    m_ionicModel->setStimulation(m_stim);  // set the stimulation
    
    //Initialize ionic model
    //________________________________
    
    m_list_ODE = m_ionicModel->getODEname();
    m_list_Current = m_ionicModel->getCurrentname();
    
    
    m_vecSeqODE.resize(m_list_ODE.size());
    m_vecSeqCurrent.resize(m_list_Current.size());
    
    
    for(unsigned int ie=0 ; ie<m_list_ODE.size() ; ie++){
      VecCreateSeq(PETSC_COMM_SELF, p_dof->numDof()/m_numComp, &m_vecSeqODE[ie]);
    }
    for(unsigned int ic=0 ; ic<m_list_Current.size() ; ic++){
      VecCreateSeq(PETSC_COMM_SELF, p_dof->numDof()/m_numComp, &m_vecSeqCurrent[ic]);
    }
    
    //________________________________
    //________________________________

    VecDestroy(&Vinit);
    
    // ionic current
    VecCreate(PETSC_COMM_WORLD,&m_ion);
    VecSetSizes(m_ion,m_localsize,p_dof->numComp()*p_mesh->nv);
    VecSetFromOptions(m_ion);
    
    
    if(m_writeODE){
      m_vecODE = m_ionicModel->getODE();
      gatherODE();
      if(m_rankProc==0){
        InitEnsight(m_name,"ODE",p_mesh,m_list_ODE,m_Tmax,m_dt,m_Nmax,m_freq_writeSolution);
        writeEnsight("ODE",0,0.,m_vecSeqODE,m_list_ODE,p_dof);
      }
    }
    if(m_writeCurrent){
      m_vecCurrent = m_ionicModel->getCurrent();
      gatherCurrent();
      if(m_rankProc==0){
        InitEnsight(m_name,"Current",p_mesh,m_list_Current,m_Tmax,m_dt,m_Nmax,m_freq_writeSolution);
        writeEnsight("Current",0,0.,m_vecSeqCurrent,m_list_Current,p_dof);
      }
    }
    
  }
  
#pragma endregion
  
  //########################//
  //########################//
  
  
  //###############################//
  // New value of Vm in ionicModel //
  //###############################//
  
#pragma region SET NEW VM VALUE IN IONIC MODEL
  
  void Bidomainproblem::setVm_new(){
    
    Vec Vm;
    VecCreate(PETSC_COMM_WORLD,&Vm);
    VecSetSizes(Vm,m_localsize/m_numComp,p_mesh->nv);
    VecSetFromOptions(Vm);

    allgather();
    
    for(int i=0 ; i<m_localsize/m_numComp ; i++){
      
      double value;
      platformInt idofloc = m_minOwnership+i*m_numComp;
      AOApplicationToPetsc(p_dof->ao(), 1, &idofloc);
      VecGetValues(m_vecSeqTmp,1,&idofloc,&value);  // m_vecSeqTmp created in the class Problem
      
      platformInt pos = i+m_minOwnership/m_numComp;
      double VmForIM = 0.;
      
      if(m_parameters.splitting){
        VmForIM = value;
      }
      else{
        if(m_ionicmodel=="FHN" || m_ionicmodel=="MV" || m_ionicmodel=="PCBE"){
          VmForIM = convertVmForIonicModel(value);
        }
        else{
          VmForIM = value;
        }
      }
      
      VecSetValues(Vm,1,&pos,&VmForIM, INSERT_VALUES);

    }
    
    VecAssemblyBegin(Vm);
    VecAssemblyEnd(Vm);
    
    m_ionicModel->setVm(Vm);
    
    VecDestroy(&Vm);
    
  }
  
#pragma endregion
  
  //###############################//
  //###############################//

#pragma region RESCALE VM (PHENOMENOLOGICAL IONIC MODEL)
  
  double Bidomainproblem::convertVmForIonicModel(double VmSol){

    double Vscaled = (VmSol-m_vmin)/m_dv;
    return Vscaled;

  }
  
  double Bidomainproblem::convertVmForSolution(double VmScaled){

    double VSol = (VmScaled*m_dv)+m_vmin;
    return VSol;

  }
  
#pragma endregion
  //###############################//
  //###############################//
  
#pragma region VM Heterogeneity initial condition
  
  int Bidomainproblem::getInitialCondition(const double cell){
    if (cell<1./3.){
      return 0;
    }
    else if (cell>1./3. and cell<2./3.){
      return 1;
    }
    else {
      return 2;
    }
  }


  void Bidomainproblem::readHeterogeneityField(vector<double> &heteroVec, string pID){
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
    //cout << "heteroVec.size() = " <<  heteroVec.size() << " - should be 3845" << endl;
  }
  
  
#pragma endregion
  
  //#######################################################################################//
  // Compute ODE (+Vm if splitting), Current and Source Term (Electrode or Electrode+Iapp) //
  //#######################################################################################//
  
#pragma region SOLVE ODE (IONIC MODEL)
  
  Vec Bidomainproblem::computeODEpart(double ttime, Vec& V_split, bool actualize){   // this function is called by the function forward (in problem.cpp)
    
    if(m_parameters.splitting){
      if(actualize){
        setVm_new();
        m_ionicModel->actualizeCVODE(ttime-(1.-m_theta)*m_dt);
      }
    }
    else{
      setVm_new();
    }
    
    m_ionicModel->solveODE(ttime);
    if(!m_parameters.splitting){
      m_ionicModel->compute_I_ion();  // computing the ionic current (call the function implemented inside ionicModel.hpp)
    }
    
    m_vecODE = m_ionicModel->getODE();
    //check_variables(m_vecODE);
    m_vecCurrent = m_ionicModel->getCurrent();
    
    if(!m_parameters.splitting){
      output();
    }
    else{
      if(actualize){
        output();
      }
    }
    
    // Electrodes
    std::vector<double> Ielsource; // I_electrodes

    if(m_parameters.splitting){
      electrodeSimulation(actualize,Ielsource); // function defined here (below) to compute Ielsource=I_el/area_el
    }
    else{
      electrodeSimulation(true,Ielsource);
    }
    //_________________________________________

    VecZeroEntries(m_ion);
    

    for(int i=0 ; i<m_localsize/m_numComp ; i++){
      
      platformInt posglob = (m_minOwnership/m_numComp)+i;
      double value;

      if(m_parameters.splitting){
        VecGetValues(m_ionicModel->getODE().back(),1,&posglob,&value);
      }
      else{
        VecGetValues(m_ionicModel->getI_ion(),1,&posglob,&value);
      }
      
      platformInt pos = (m_minOwnership/m_numComp)+i;
      platformInt idof = pos*p_dof->numComp();
      AOApplicationToPetsc(p_dof->ao(), 1, &idof);
      
      if( (m_ionicmodel=="FHN" || m_ionicmodel=="MV" || m_ionicmodel=="PCBE") && m_parameters.splitting==false){
        value = value*m_dv*m_Cm;
      }
      
      if(m_parameters.splitting){
        VecSetValues(V_split,1,&idof,&value, INSERT_VALUES);
      }
      else{
        VecSetValues(m_ion,1,&idof,&value, INSERT_VALUES);
      }
      
    }
    
    VecAssemblyBegin(m_ion);
    VecAssemblyEnd(m_ion);

    if(m_ionicmodel != "Paci"){//TO be in the same condition for Paci (for validation)

      for(int i=0 ; i<m_localsize/m_numComp ; i++){

        platformInt posglob = (m_minOwnership/m_numComp)+i;
        platformInt pos = (m_minOwnership/m_numComp)+i;
        platformInt idof = pos*p_dof->numComp();// Vm
        
        AOApplicationToPetsc(p_dof->ao(), 1, &idof);
        
        const double iapp = m_stim.value_Iapp(posglob,m_time);
        //VecSetValues(m_ion,1,&idof,&iapp, ADD_VALUES);
        
      }
      
    }
    
    
    if(!m_parameters.splitting){
      VecAssemblyBegin(m_ion);
      VecAssemblyEnd(m_ion);
    }
    
    // Neumann Conditions (electrodes)
    electrodeNeumannConditions(Ielsource);
    //____________________________________________
    
    if(m_parameters.splitting){
      VecAssemblyBegin(V_split);
      VecAssemblyEnd(V_split);
    }
    
    return m_ion;
  }
  
#pragma endregion
  
  //#######################################################################################//
  //#######################################################################################//
  
  
  //##################################//
  // Call functions for posttreatment //
  //##################################//
  
#pragma region TO WRITE SOLUTION
  
  void Bidomainproblem::output(){
    
    if(m_writeODE && m_numIt%m_freq_writeSolution==0){
      gatherODE();
      const int num = m_numIt/m_freq_writeSolution;
      if(m_rankProc==0){
        writeEnsight("ODE",num,m_time,m_vecSeqODE,m_list_ODE,p_dof);
      }
    }
    if(m_writeCurrent && m_numIt%m_freq_writeSolution==0){
      gatherCurrent();
      const int num = m_numIt/m_freq_writeSolution;
      if(m_rankProc==0){
        writeEnsight("Current",num,m_time,m_vecSeqCurrent,m_list_Current,p_dof);
      }
    }
    
  }
  
#pragma endregion
  
  //##################################//
  //##################################//
  
  
  //############################//
  // GROUNDS BOUNDARY CONDITION //
  //############################//

#pragma region GROUND MEA MODEL
  void Bidomainproblem::groundBC(){

    const double sigma_val = m_sigmai + m_sigmae;
    const double factor = 1.0;
    
    for(map< int, std::vector<int> >::const_iterator it=p_mesh->vertexPerLabel.begin() ; it!=p_mesh->vertexPerLabel.end() ; it++){
      int label = it->first;
      std::vector<int> noeuds = it->second;
      
      if(label==5){//Ce label doit correspondre au label des masses
        for(unsigned int in=0 ; in<noeuds.size() ; in++){
          
          platformInt idofUe = noeuds[in]*p_dof->numComp()+1; // +1 stands for the Ue component
          AOApplicationToPetsc(p_dof->ao(), 1, &idofUe);
          
          double groundVal = 0.;
          VecGetValues(m_vecSeqTmp,1,&idofUe,&groundVal);//Get Ue
          groundVal = factor*sigma_val;

          VecSetValues(m_ion,1,&idofUe,&groundVal, ADD_VALUES);
          
        }
      }
    }
    
    VecAssemblyBegin(m_ion);
    VecAssemblyEnd(m_ion);
    
  }
#pragma endregion

  //############################//
  //############################//
  
  //#############################################//
  // Methods related to the electrode simulation //
  //#############################################//
  
#pragma region ELECTRODE SIMULATION
  
  void Bidomainproblem::electrodeSimulation(bool actualize, std::vector<double>& Ielsource){
    
    if(actualize){
      computeElectrodesMean(); // mean over the electrodes of Vm and ue: m_Vmmean and m_Uemean
      //writeElectrodesSolution();
      m_electrodes.set_Umean(m_Uemean); //setting the mean over the electrodes
      m_electrodes.forwardElectrodes(); //solve the ODE for each electrode, computing Iel
      
      std::vector<double> Umes;
      m_electrodes.getUmes(Umes);
      
      //writeUemes(Umes);
      for(int ie=0 ; ie<m_nbEl ; ie++){
        m_UmesVec[ie].push_back(Umes[ie]);
        m_VmVec[ie].push_back(m_Vmmean[ie]);
      }
      if(m_parameters.ionicmodel == "OHaraRudy"){
        m_CaiVec.push_back(m_CaimeanWell);
      }
    }
    //To have source term for the electrode simulation
    m_electrodes.computeSourceIel(Ielsource);   // function from electrodes.cpp (compute I_el/area_el)
  }
  
  void Bidomainproblem::electrodeNeumannConditions(const std::vector<double>& Ielsource){
    
    for(map< int, std::vector<int> >::const_iterator it=p_mesh->vertexPerLabel.begin() ; it!=p_mesh->vertexPerLabel.end() ; it++){
      int label = it->first;
      std::vector<int> noeuds = it->second;
      
      if(label>=10){
        for(unsigned int in=0 ; in<noeuds.size() ; in++){
          
          platformInt idofUe = noeuds[in]*p_dof->numComp()+1; // +1 stands for the Ue component
          AOApplicationToPetsc(p_dof->ao(), 1, &idofUe);
          
          const double Iel_Neumann = Ielsource[label-10]/m_sigmae;
          VecSetValues(m_ion,1,&idofUe,&Iel_Neumann, ADD_VALUES);
          
        }
      }
    }
    
    VecAssemblyBegin(m_ion);
    VecAssemblyEnd(m_ion);
    
  }
  
#pragma endregion
  
  //#############################################//
  //#############################################//
  
  
  //###########################################//
  // Write the mean solution of each electrode //
  //###########################################//
  
#pragma region WRITE VM UE MEAN ON ELECTRODE

  void Bidomainproblem::writeElectrodesSolution_headers(){
    
    if(m_rankProc==0){
      
      vector<string> dummy;
      dummy.push_back("value");
      for(int ie=0 ; ie<m_nbEl ; ie++){
        
        std::stringstream index;
        index << ie;
        
        string filenameUe = "./Ue/Ue" + index.str() + "_p" + m_procId + ".txt";
        string filenameVm = "./Vm/Vm" + index.str() + "_p" + m_procId + ".txt";
        string filenameUmes = "./Umes/Umes" + index.str() + "_p" + m_procId + ".txt";
        
        writeSolution_0d_Headers(dummy,filenameUe,m_parameters.binaryPostProc);
        writeSolution_0d_Headers(dummy,filenameVm,m_parameters.binaryPostProc);
        writeSolution_0d_Headers(dummy,filenameUmes,m_parameters.binaryPostProc);
        

      }

      if(m_parameters.ionicmodel == "OHaraRudy"){
        string filenameCai = "./Cai/CaiWell_p" + m_procId + ".txt";
        writeSolution_0d_Headers(dummy,filenameCai,m_parameters.binaryPostProc);
      }
    }
    
  }

  void Bidomainproblem::writeElectrodesSolution(){
    
    if(m_rankProc==0){
      
      for(int ie=0 ; ie<m_nbEl ; ie++){
        
        std::stringstream index;
        index << ie;
        
        string filenameUe = "./Ue/Ue" + index.str() + "_p" + m_procId + ".txt";
        string filenameVm = "./Vm/Vm" + index.str() + "_p" + m_procId + ".txt";
        
        writeSolution_0d(m_time,m_Uemean[ie],filenameUe);
        writeSolution_0d(m_time,m_Vmmean[ie],filenameVm);

      }
      
      for(int ie=0 ; ie<m_nbEl ; ie++){
        SolutionElectrodes[ie].push_back(m_Uemean[ie]);
        SolutionVm[ie].push_back(m_Vmmean[ie]);
      }
      
      if(m_parameters.ionicmodel == "OHaraRudy"){
        string filenameCaiWell = "./Cai/CaiWell_p" + m_procId + ".txt";
        writeSolution_0d(m_time, m_CaimeanWell, filenameCaiWell);
      }
    }
    
  }
  
#pragma endregion
  
  //###########################################//
  //###########################################//


  //############//
  // Write Umes //
  //############//
  
#pragma region WRITE MEASURED POTENTIAL ON ELECTRODE
  
  void Bidomainproblem::writeUemes(const std::vector<double>& Umes){
    
    if(m_rankProc==0){
      
      for(int ie=0 ; ie<m_nbEl ; ie++){
        
        std::stringstream index;
        index << ie;
        
        string filenameUmes = "./Umes/Umes" + index.str() + "_p" + m_procId + ".txt";
        writeSolution_0d(m_time,Umes[ie],filenameUmes);
        
      }
    }

  }

  void Bidomainproblem::writeUemesAtTheEnd(){
    
    if(m_rankProc==0){
      
      for(int ie=0 ; ie<m_nbEl ; ie++){
        stringstream fileName;
        fileName << "./Umes/Umes" << ie << "_p" << m_procId << ".txt";
        ofstream outFile(fileName.str().c_str());
        outFile<<scientific;
        outFile.precision(15);
        vector<double> umes = m_UmesVec[ie];
        size_t n = umes.size(); 
        for (int i=0; i<n; i++){
          outFile << umes[i] << " ";
        }
        outFile << "\n";
        outFile.close();


        stringstream fileName2;
        fileName2 << "./Vm/Vm" << ie << "_p" << m_procId << ".txt";
        ofstream outFile2(fileName2.str().c_str());
        outFile2<<scientific;
        outFile2.precision(15);
        vector<double> vm = m_VmVec[ie];
        n = vm.size(); 
        for (int i=0; i<n; i++){
          outFile2 << vm[i] << " ";
        }
        outFile2 << "\n";
        outFile2.close();
        
      }
      if(m_parameters.ionicmodel == "OHaraRudy"){
        stringstream fileNameCai;
        fileNameCai << "./Cai/CaiWell_p" << m_procId << ".txt";
        ofstream outFile(fileNameCai.str().c_str());
        outFile<<scientific;
        outFile.precision(15);
        size_t n = m_CaiVec.size(); 
        for (int i=0; i<n; i++){
          outFile << m_CaiVec[i] << " ";
        }
        outFile << "\n";
        outFile.close();
      }
      
    }

  }
  
#pragma endregion
  
  //############//
  //############//
  
  
  //#########################################//
  // Compute Mean quantity on each electrode //
  //#########################################//
  
#pragma region COMPUTE MEAN VALUES ON ELECTRODE
  
  void Bidomainproblem::computeElectrodesMean(){
    
    if(m_time==m_dt){

      m_labels.clear();
      m_measure.clear();
      
      for(map< int, std::vector<int> >::const_iterator it=p_mesh->vertexPerLabel.begin() ; it!=p_mesh->vertexPerLabel.end() ; it++){
        int label = it->first;
        
        if(label>=10){//We are on an electrode
          m_labels.push_back(label);
        }
      }
      
      m_measure.resize(m_labels.size());
      for(int ie=0 ; ie<p_mesh->nt ; ie++){
        for(unsigned int ilab=0 ; ilab<m_labels.size() ; ilab++){
          if(p_mesh->triangles[ie].ref == m_labels[ilab]){
            m_measure[m_labels[ilab]-10] += p_mesh->triangles[ie].area;
          }
        }
      }
    }
    
    m_Uemean.clear();
    m_Vmmean.clear();
    
    m_Uemean.resize(m_labels.size());
    m_Vmmean.resize(m_labels.size());

    for(map< int, std::vector<int> >::const_iterator it=p_mesh->vertexPerLabel.begin() ; it!=p_mesh->vertexPerLabel.end() ; it++){
      
      int label = it->first;
      
      if(label>=10){//We are on an electrode
        
        std::vector<int> noeuds = it->second;
        double sommeUe = 0., sommeVm = 0.;
        
        for(unsigned int in=0 ; in<noeuds.size() ; in++){
          
          platformInt idofVm = noeuds[in]*p_dof->numComp();
          platformInt idofUe = noeuds[in]*p_dof->numComp()+1;
          
          AOApplicationToPetsc(p_dof->ao(), 1, &idofVm);
          AOApplicationToPetsc(p_dof->ao(), 1, &idofUe);
          
          double Vmnoeud,Uenoeud;
          VecGetValues(m_vecSeqTmp,1,&idofVm,&Vmnoeud);//Sol
          VecGetValues(m_vecSeqTmp,1,&idofUe,&Uenoeud);//Sol
          
          sommeUe += Uenoeud;
          sommeVm += Vmnoeud;
          
        }

        m_Uemean[label-10] = sommeUe/noeuds.size();
        m_Vmmean[label-10] = sommeVm/noeuds.size();
        
      }
    }
    
    //Mean on the whole well
    double sumCaiWell = 0.;
    int nbNoeuds = 0;
    if(m_parameters.ionicmodel == "OHaraRudy"){
      
      // for(map< int, std::vector<int> >::const_iterator it=p_mesh->vertexPerLabel.begin() ; it!=p_mesh->vertexPerLabel.end() ; it++){
      //   int label = it->first;
      //   std::vector<int> noeuds = it->second;
        
      //   for(unsigned int in=0 ; in<noeuds.size() ; in++){
          
      //     platformInt idofCai = noeuds[in]*p_dof->numComp();
      //     AOApplicationToPetsc(p_dof->ao(), 1, &idofCai);
          
      //     double Cainoeud;
      //     VecGetValues(m_vecODE[36],1,&idofCai,&Cainoeud);
      //     //VecGetValues(m_vecODE[38],1,&idofCai,&Cainoeud);
      //     sumCaiWell += Cainoeud;
          
      //   }
      //   nbNoeuds += noeuds.size();
      // }
      //cout << m_time << endl;
      if(m_time<=m_dt){
        //m_CaimeanWell = m_variables.m_InitialCond[getInitialCondition(0)][36];
        m_CaimeanWell = m_variables.Cai;
      }
      else{
        m_CaimeanWell = computeCaiMean(); //sumCaiWell/double(nbNoeuds);
      }
    }
    
  }
  
  
  
  double Bidomainproblem::computeCaiMean(){
    double val;
    PetscInt dummy;
    Vec vecCaiSol;
    gatherODE();
    //Vec vecCai = m_vecODE[36];//36 for ORd model
    //VecGetSize(vecCai, &dummy);
    //cout << "vecCai size = " << dummy << endl;
    //VecGetSize(p_vec[0], &dummy);
    //cout << "vecCai size = " << dummy << endl;
    //cout << "avant vecdupltire" << dummy << endl;
    VecDuplicate(p_vec[0], &vecCaiSol);
    //cout << "after vecduplicate" << dummy << endl;
    VecSet(vecCaiSol, 0.);
    //VecGetSize(vecCaiSol, &dummy);
    //cout << "vecCaiSol size = " << dummy << endl;

    for(map< int, std::vector<int> >::const_iterator it=p_mesh->vertexPerLabel.begin() ; it!=p_mesh->vertexPerLabel.end() ; it++){
      std::vector<int> noeuds = it->second;
      for(unsigned int in=0 ; in<noeuds.size() ; in++){
          
        platformInt idofVm = noeuds[in]*p_dof->numComp();
        platformInt idofUe = noeuds[in]*p_dof->numComp(); //+1;
        AOApplicationToPetsc(p_dof->ao(), 1, &idofVm);
        AOApplicationToPetsc(p_dof->ao(), 1, &idofUe);
        //cout << idofVm << " " << idofUe << endl;
        idofVm /=2;
        VecGetValues(m_vecODE[36],1,&idofVm,&val);
        VecSetValues(vecCaiSol,1,&idofUe,&val,INSERT_VALUES);
      }
      VecAssemblyBegin(vecCaiSol);
      VecAssemblyEnd(vecCaiSol);
    }

    VecGetSize(vecCaiSol, &dummy);
    //cout << "vecCaiSol size = " << dummy << endl;
    
    PetscScalar integralOfCai, volume, averaged;
    Vec vecOfOnes, tmpVec;
    VecDuplicate(vecCaiSol, &tmpVec);
    VecDuplicate(vecCaiSol, &vecOfOnes);
    VecSet(vecOfOnes, 1.0);
    
    MatMult(p_mat[1], vecCaiSol, tmpVec);
    VecDot(vecOfOnes, tmpVec, &integralOfCai);

    //Compute domain volume
    MatMult(p_mat[1], vecOfOnes, tmpVec);
    VecDot(vecOfOnes, tmpVec, &volume);
    averaged = integralOfCai/volume;
    //cout << "mean Cai = " << averaged << endl;
    //cout << "before suspect" << endl;
    VecDestroy(&vecCaiSol);
    VecDestroy(&vecOfOnes); VecDestroy(&tmpVec);
    //cout << "after suspect" << endl;
    return averaged; 
  }
  
  
#pragma endregion
  
  //#########################################//
  //#########################################//
  
#pragma region GATHER (ODE-CURRENT)
  
  void Bidomainproblem::gatherODE(){
    
    for(unsigned int ie=0 ; ie<m_list_ODE.size() ; ie++){
      Vec vout;
      VecScatter vecscat;
      
      VecScatterCreateToAll(m_vecODE[ie],&vecscat,&vout);
      VecScatterBegin(vecscat,m_vecODE[ie],vout,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(vecscat,m_vecODE[ie],vout,INSERT_VALUES,SCATTER_FORWARD);
      
      VecScatterDestroy(&vecscat);
      VecCopy(vout,m_vecSeqODE[ie]);
      VecDestroy(&vout);
    }
    
  }
  
  void Bidomainproblem::gatherCurrent(){
    
    for(unsigned int ic=0 ; ic<m_list_Current.size() ; ic++){
      Vec vout;
      VecScatter vecscat;
      
      VecScatterCreateToAll(m_vecCurrent[ic],&vecscat,&vout);
      VecScatterBegin(vecscat,m_vecCurrent[ic],vout,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(vecscat,m_vecCurrent[ic],vout,INSERT_VALUES,SCATTER_FORWARD);
      
      VecScatterDestroy(&vecscat);
      VecCopy(vout,m_vecSeqCurrent[ic]);
      VecDestroy(&vout);
    }

  }
  
#pragma endregion
  
  //#################//
  //#################//
  
#pragma region INITIALIZE ELECTRODE
  
  void Bidomainproblem::initializeElectrodeParameters(int nbEl, const vector<double>& Ri, const vector<double>& Rel, const vector<double>& Cel, const vector<double>& Diameters){
    
    m_Ri.clear();
    m_Rel.clear();
    m_Cel.clear();
    m_diametres.clear();
    
    if(Ri.size()==nbEl){
      m_Ri = Ri;
    }
    else{
      for(int ie=0 ; ie<nbEl ; ie++){
        m_Ri.push_back(Ri[0]);
      }
    }

    if(Rel.size()==nbEl){
      m_Rel = Rel;
    }
    else{
      for(int ie=0 ; ie<nbEl ; ie++){
        m_Rel.push_back(Rel[0]);
      }
    }
    
    if(Cel.size()==nbEl){
      m_Cel = Cel;
    }
    else{
      for(int ie=0 ; ie<nbEl ; ie++){
        m_Cel.push_back(Cel[0]);
      }
    }

    if(Diameters.size()==nbEl){
      m_diametres = Diameters;
    }
    else{
      for(int ie=0 ; ie<nbEl ; ie++){
        m_diametres.push_back(Diameters[0]);
      }
    }

  }

#pragma endregion

	//#################//
	//#################//

}
