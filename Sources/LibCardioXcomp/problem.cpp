#include "stdafx.h"

#include "problem.hpp"

#include <ctime>

using namespace std;

namespace cardioxcomp {
  
  //#############//
  // CONSTRUCTOR //
  //#############//
  
#pragma region VOID CONSTRUCTOR
  
  Problem::Problem(){
  }
  
#pragma endregion
  
#pragma region BASIC CONSTRUCTOR
  
  Problem::Problem(string meshfile,int numComp,int numMat,int numVec,int numProc,int rankProc, bool isHetero, string pythonProc, bool cellHetero):
    m_numComp(numComp),m_numMat(numMat),m_numVec(numVec),m_numProc(numProc),m_rankProc(rankProc),m_isHetero(isHetero),m_pythonProc(pythonProc),m_cellHetero(cellHetero)
  {
    p_mesh = new Mesh(meshfile,m_rankProc);
    p_dof = new Dof(*p_mesh,m_numComp);
    p_dof->initializePattern(numProc,rankProc);
    p_dof->partitionDof(numProc,rankProc,PETSC_COMM_WORLD);
    p_mat = new Mat[m_numMat];
    p_vec = new Vec[m_numVec]; //table of Vec: first is the solution, second is the rhs integrated, third is the rhs function
    p_matElem = new MatElem[m_numMat]; //table of matrices: first is the bidomain problem matrix, second is the mass matrix
    p_vecElem = new VecElem[m_numVec];
    for(int i=0;i<m_numMat;i++){
      p_matElem[i].init(m_numComp);
      p_matElem[i].zero();
    }
    for(int i=0;i<m_numVec;i++){
      p_vecElem[i].init(m_numComp);
      p_vecElem[i].zero();
    }
    p_localMesh = new Mesh(*p_mesh,p_dof->eltPart(),m_rankProc,m_loc2GlobElem);

#ifdef PETSC_3_6
    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, m_loc2GlobElem.size(),1, &m_loc2GlobElem[0], PETSC_COPY_VALUES, &m_mappingElem);
#else
    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, m_loc2GlobElem.size(), m_loc2GlobElem.data(), PETSC_COPY_VALUES, &m_mappingElem);
#endif

    m_initMappingElem=true;
    p_dof->allocateMatrixVec(m_numProc, m_rankProc, p_mat, m_numMat, p_vec, m_numVec);
    VecCreateSeq(PETSC_COMM_SELF, p_dof->numDof(), &m_vecSeqTmp);
    
    platformInt dummy;
    VecGetLocalSize(p_vec[0],&m_localsize);
    VecGetOwnershipRange(p_vec[0],&m_minOwnership,&dummy);
    
    VecCreate(PETSC_COMM_WORLD,&m_extrapolate);
    VecSetSizes(m_extrapolate,m_localsize,m_numComp*p_mesh->nv);
    VecSetFromOptions(m_extrapolate);
    
    VecSet(m_extrapolate,0.);
    
    VecAssemblyBegin(m_extrapolate);
    VecAssemblyEnd(m_extrapolate);
    
    m_time = 0.;
    m_numIt = 0;

    if (m_cellHetero){
      readInitComp();
    }
    
  }
  
#pragma endregion
  
  //#############//
  //#############//
  
  
  //############//
  // DESTRUCTOR //
  //############//
  
#pragma region DESTRUCTOR

#ifdef _WIN32
  __declspec(dllexport) 
#endif
  Problem::~Problem(){
    
    if(m_initMappingElem) ISLocalToGlobalMappingDestroy(&m_mappingElem);
    for(int i=0;i<m_numMat;i++) MatDestroy(&p_mat[i]);
    for(int i=0;i<m_numVec;i++) VecDestroy(&p_vec[i]);
    VecDestroy(&m_vecSeqTmp);
    
    VecDestroy(&m_source);
    KSPDestroy(&m_ksp);
    //VecDestroy(&m_extrapolate);
    
    delete p_dof;
    delete p_matElem;
    delete p_vecElem;
    
    delete p_mat;
    delete p_vec;
    
  }
  
#pragma endregion
  
  //############//
  //############//
  
  
  //#################//
  // FINALIZE MATVEC //
  //#################//
  
#pragma region FINALIZE MAT-VEC
  
  void Problem::finalizeMatVec(){
    for(int i=0;i<m_numMat;i++){
      MatAssemblyBegin(p_mat[i],MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(p_mat[i],MAT_FINAL_ASSEMBLY);
    }
    
    for(int i=0;i<m_numVec;i++){
      VecAssemblyBegin(p_vec[i]);
      VecAssemblyEnd(p_vec[i]);
    }
  }
  
#pragma endregion
  
  //#################//
  //#################//
  
  
  //########//
  // SOLVER //
  //########//
  
#pragma region BUILD SOLVER
  
  void Problem::buildSolver(){
    
    KSPCreate(PETSC_COMM_WORLD, &m_ksp);
    
    KSPSetType(m_ksp,m_solver.c_str());
    KSPGMRESSetRestart(m_ksp,200);
    
    KSPGetPC(m_ksp,&m_pc);
    PCSetFromOptions(m_pc);
    PCSetType(m_pc,m_preconditioner.c_str());
    
#ifndef PETSC_3_6
    if(m_preconditioner == "lu"){
      PCFactorSetMatSolverPackage(m_pc,MATSOLVERSUPERLU_DIST);
    }
#endif
    
    KSPSetTolerances(m_ksp,m_edpReltol,m_edpAbstol,PETSC_DEFAULT,m_edpMaxiter);
    
    if(m_preconditionerOption == "SAME_NONZERO_PATTERN"){
      m_setPreconditionerOption = SAME_NONZERO_PATTERN;
    }
    else if(m_preconditionerOption == "SAME_PRECONDITIONER"){
#ifndef PETSC_3_6
      m_setPreconditionerOption = SAME_PRECONDITIONER;
#endif
    }
    else if(m_preconditionerOption == "DIFFERENT_NONZERO_PATTERN"){
      m_setPreconditionerOption = DIFFERENT_NONZERO_PATTERN;
    } else {
      std::cout << "Error in config.lua file (preconditionerOptions): SAME_NONZERO_PATTERN, SAME_PRECONDITIONER or DIFFERENT_NONZERO_PATTERN" << std::endl;
      exit(0);
    }
    
    
    
  }
  
#pragma endregion
  
  //########//
  //########//
  
  
  //###############//
  //  Remove Mean  //
  //###############//
  
#pragma region REMOVE MEAN
  
  void Problem::removeMean(int comp)
  {
    allgather();
    int size = p_dof->numDof()/p_dof->numComp();
    //int size = p_mesh->nv;
     
    double Ue_mean;
    double sumUe = 0.;
    double minVm = 10000;
    double maxVm = -10000;
    for (int i = 0; i < size; i++){
      
      platformInt idof = i*p_dof->numComp()+comp;//comp = 0 -> idofVm, comp = 1 -> idofUe
      platformInt idofVm = i*p_dof->numComp();
      AOApplicationToPetsc(p_dof->ao(), 1, &idof);
      AOApplicationToPetsc(p_dof->ao(), 1, &idofVm);
      
      double Ue_node;
      double Vm_node;
      VecGetValues(m_vecSeqTmp,1,&idof,&Ue_node);// get values of Vm if comp = 0, Ue if comp = 1
      VecGetValues(m_vecSeqTmp,1,&idofVm,&Vm_node);
      sumUe += Ue_node;
      if(Vm_node > maxVm)
        maxVm = Vm_node;
      if(Vm_node<minVm)
        minVm = Vm_node;
    }
    Ue_mean = sumUe/size;// mean value on all the MEA
    if (m_rankProc==0)
      {
        std::cout<<" The mean value of Extracellular is "<< Ue_mean<<std::endl;
        std::cout<<" minVm, maxVm is "<< minVm<<" "<<maxVm<<std::endl;
      }

    
    for(int i=0 ; i<m_localsize/m_numComp ; i++){
      
      double Ue_node;
      platformInt idofloc = i*m_numComp+1;
      platformInt idof;

#ifdef PETSC_3_6
      ISLocalToGlobalMappingApply(p_dof->locToGlob(),1,&idofloc,&idof);
#else
      ISLocalToGlobalMappingApply(p_dof->return_mappingNodes(),1,&idofloc,&idof);
#endif

      AOApplicationToPetsc(p_dof->ao(), 1, &idof);      // ordering application

      VecGetValues(p_vec[0],1,&idof,&Ue_node);   // extrapolated potential inside the variable value
      Ue_node=Ue_node-Ue_mean;
      VecSetValues(p_vec[0],1,&idof,&Ue_node, INSERT_VALUES); // potential inside Vinit 
    }
  
    VecAssemblyBegin(p_vec[0]);
    VecAssemblyEnd(p_vec[0]);
    
  }
  
#pragma endregion

  //###############//
  //###############//
  
  
  //############//
  // SOLVE AX=B //
  //############//
  
#pragma region SOLVE SYSTEM AX=B
  
  void Problem::solve(){
    
    Mat& matrix = p_mat[0]; // this is the bidomain problem matrix (for the lhs)
    
    Mat& matrixMass = p_mat[1]; // this is the mass matrix (used in the rhs)
    
    Vec& sol = p_vec[0]; // solution vector
    Vec& rhs = p_vec[1]; // rhs vector integrated
    Vec& fct = p_vec[2]; // function of the rhs vector
    
    MatMult(matrixMass,fct,rhs); // matrix-vector product: rhs = matrixMass*fct 
    PCFactorSetReuseFill(m_pc, PETSC_TRUE);
    
    // KSP solver
#ifdef PETSC_3_6    
    KSPSetReusePreconditioner(m_ksp,PETSC_TRUE);
    KSPSetOperators(m_ksp,matrix,matrix);
#else
    KSPSetOperators(m_ksp,matrix,matrix,m_setPreconditionerOption);
#endif
    KSPSetFromOptions(m_ksp);
    KSPSolve(m_ksp,rhs,sol); // solve the system: matrix*sol=rhs obtaining sol
  }
  
#pragma endregion
  
  //############//
  //############//
  
  
  //##########//
  // POSTPROC //
  //##########//
  
#pragma region POST-PROCESSING
  
  void Problem::postprocScalar(const string& filename,int ivec,int icomp){
    Vec& sol = p_vec[ivec];
    
    Vec vout;
    VecScatter vecscat;
    VecScatterCreateToAll(sol,&vecscat,&vout);
    VecScatterBegin(vecscat,sol,vout,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(vecscat,sol,vout,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&vecscat);
    VecCopy(vout,m_vecSeqTmp);
    VecDestroy(&vout);
    
    platformInt* iout = new platformInt[p_dof->numDof()];
    double* valueVec = new double[p_dof->numDof()];
    for (int i=0; i<p_dof->numDof(); i++) iout[i] = i;
    AOApplicationToPetsc(p_dof->ao(),p_dof->numDof(),iout);
    VecGetValues(m_vecSeqTmp,p_dof->numDof(),iout,valueVec);
    delete [] iout;
    
    if(m_rankProc == 0){
      int size =p_dof->numDof()/p_dof->numComp();
      FILE* pFile;
      pFile = fopen(filename.c_str(),"w" );
      fprintf(pFile,"MeshVersionFormatted \n 2 \n");
      fprintf(pFile,"Dimension \n 2 \n");
      fprintf(pFile,"SolAtVertices \n %d \n", size);
      fprintf(pFile,"1 1\n");
      for ( int i = 0; i < size; i++)
        fprintf(pFile, "%le \n", valueVec[i*p_dof->numComp()+icomp]);
      fprintf(pFile,"End");
      fclose(pFile);
    }
  }
  
#pragma endregion
  
  //##########//
  //##########//
  
  
  //#########################//
  // FORWARD: time iteration //
  //#########################//
  
#pragma region CLASSIC FORWARD (WITHOUT SPLITTING)
  
  void Problem::forward(){
    
    
    //------------------------
    // INITIALIZATION
    //------------------------
      
    if(m_time==0.){
      
      
      m_bdfEdp.defineOrder(m_orderbdf,1);
      
      for (int i=0 ; i<m_localsize/m_numComp ; i++) {
        
        for(int ivar=0 ; ivar<m_numComp ; ivar++){
          
          platformInt pos = i+m_minOwnership/m_numComp;
          platformInt idof = pos*p_dof->numComp()+ivar;
          AOApplicationToPetsc(p_dof->ao(), 1, &idof);

          
          if (!m_cellHetero){
            VecSetValue(p_vec[0], idof, m_physicalInitialValues[ivar], INSERT_VALUES); // initial values in p_vec[0]
          }
          else{
            if (ivar==0){
              VecSetValue(p_vec[0], idof, m_comp1Init[pos], INSERT_VALUES); // initial values in p_vec[0]
            }
            else{
              VecSetValue(p_vec[0], idof, m_physicalInitialValues[ivar], INSERT_VALUES); // initial values in p_vec[0]
            }
          }
          
        }
      }
      
      
      VecAssemblyBegin(p_vec[0]);
      VecAssemblyEnd(p_vec[0]);
      
        // $$$$$$$$$$$$$$$ What should we do for intializing the rescaled potential: call new initialization  $$$$$$$$$$$$$$$//
        
      m_bdfEdp.initialize(p_vec[0]);
      m_bdfEdp.extrapolate(m_extrapolate);
      initializeProblem(m_extrapolate);//Create ionicModel
      //------------------------
      //------------------------
      
      allgather(); //creating the sequential vector for ensight files
      if(m_rankProc==0){
        InitEnsight(m_name,p_mesh,m_physicalvariables,m_Tmax,m_dt,m_Nmax,m_freq_writeSolution);
        writeEnsight(m_name,0,0.,m_vecSeqTmp,m_physicalvariables,p_dof);
      }
      
    }
    
    //------------------------
    // TIME > 0
    //------------------------
    
    else{
      m_bdfEdp.update(p_vec[0]);// update the solution according to bdf order: new p_vec[0]
      m_bdfEdp.extrapolate(m_extrapolate);// extrapolation of m_extrapolate with the solution
    }
    
    
    //------------------------
    // SOLVE THE PROBLEM
    //------------------------
    m_bdfEdp.computeRHSTime(m_dt); // compute the rhs at the current time
    
    if (m_time==0.){
      computeMatVec(); // compute the bidomain matrix and the mass matrix
      dirichletBC();
      buildSolver();
    }
    
    m_time += m_dt; // updating time and iterations
    m_numIt++;
    
    p_vec[0] = m_extrapolate;
    //VecCopy(m_extrapolate,p_vec[0]);
    allgather(); // creates m_vecSeqTmp with the solution (for ensight files)
    
    VecAssemblyBegin(p_vec[2]);
    VecAssemblyEnd(p_vec[2]);
    
    VecZeroEntries(p_vec[2]); // recomputing the rhs function at the current time
    VecAXPY(p_vec[2],m_coefftime,m_bdfEdp.rhs());// dtu-lapl(u) = p_vec2----> pvec_2 += rhs

    Vec Vdummy;//usefull only with splitting
    bool Vtdummy = true;//usefull only with splitting
    
    m_source = computeODEpart(m_time,Vdummy,Vtdummy); // compute the source: Iion, Iapp, Iel
    
    VecAXPY(p_vec[2],m_coeffsource,m_source); // dtu-lapl(u) = p_vec2----> pvec_2 += source
    
    finalizeMatVec();  // assembling p_vec and p_mat
    
    solve();  // compute sol
    
    if(m_rankProc==0){
      if(m_numIt%m_freq_writeSolution == 0){
        const int num = m_numIt/m_freq_writeSolution;
        writeEnsight(m_name,num,m_time-m_dt,m_vecSeqTmp,m_physicalvariables,p_dof);
      }
    }
    
  }
  
#pragma endregion
  
  //#########################//
  //#########################//
  
  
  //##############################//
  // FORWARD FOR SPLITTING METHOD //
  //##############################//
  
#pragma region SPLITTING FORWARD
  
  void Problem::forward_split(){
    /*
     Strang splitting (theta = 1/2)
     
     "An operator splitting method for solving the bidomain equations coupled to a volume conductor model for the torso"
     Joakim Sundnes, Glenn Terje Lines, Aslak Tveito
     Mathematical Biosciences 194 (2005) 233–248
    */
     
    if(m_time==0.){
      
      
      //------------------------
      // INITIALIZATION
      //------------------------
      
      m_bdfEdp.defineOrder(m_orderbdf,1);
      
      for (int i=0 ; i<m_localsize/m_numComp ; i++) {
        
        for(int ivar=0 ; ivar<m_numComp ; ivar++){
          
          platformInt pos = i+m_minOwnership/m_numComp;
          platformInt idof = pos*p_dof->numComp()+ivar;
          AOApplicationToPetsc(p_dof->ao(), 1, &idof);
          VecSetValue(p_vec[0], idof, m_physicalInitialValues[ivar], INSERT_VALUES);
        }
        
      }
      
      
      VecAssemblyBegin(p_vec[0]);
      VecAssemblyEnd(p_vec[0]);
      
      // m_bdfEdp.initialize(p_vec[0]);
       // initializeProblem(p_vec[0]);
      // Changed because we need to get the heterogeneous initial condition before initializing BDF
      initializeProblem(p_vec[0]);
      m_bdfEdp.initialize(p_vec[0]);
     

      //------------------------
      //------------------------
      
      allgather();
      if(m_rankProc==0){
        InitEnsight(m_name,p_mesh,m_physicalvariables,m_Tmax,m_dt,m_Nmax,m_freq_writeSolution);
        writeEnsight(m_name,0,0.,m_vecSeqTmp,m_physicalvariables,p_dof);
      }
 
    }
    
    if (m_time==0.){
      computeMatVec();
      dirichletBC();
      buildSolver();
      finalizeMatVec();
    }
    
    
    //-----------
    // Step 1
    //-----------
    
    m_source = computeODEpart(m_time+m_theta*m_dt,p_vec[0],false);
    
    //-----------
    //-----------
    
    m_time += m_dt;
    m_numIt++;
    
    m_bdfEdp.update(p_vec[0]);
    m_bdfEdp.computeRHSTime(m_dt);
    
    
    VecAssemblyBegin(p_vec[2]);
    VecAssemblyEnd(p_vec[2]);
    
    VecZeroEntries(p_vec[2]);
    VecAXPY(p_vec[2],m_coefftime,m_bdfEdp.rhs());// dtu-lapl(u) = p_vec2----> pvec_2 += rhs
    VecAXPY(p_vec[2],m_coeffsource,m_source);// dtu-lapl(u) = p_vec2----> pvec_2 += source (electrode modelisation)
    
    finalizeMatVec();
    
    //-----------
    // Step 2
    //-----------
    
    solve();
    
    //-----------
    //-----------
    

    //-----------
    // Step 3
    //-----------
    
    m_source = computeODEpart(m_time,p_vec[0],true);
    
    
    allgather();
    if(m_rankProc==0){
      if(m_numIt%m_freq_writeSolution == 0){
        const int num = m_numIt/m_freq_writeSolution;
        writeEnsight(m_name,num,m_time,m_vecSeqTmp,m_physicalvariables,p_dof);
      }
    }
    
    //-----------
    //-----------
    
  }
  
#pragma endregion
  
  //##############################//
  //##############################//
  
  
  //###########################//
  //###########################//
  
#pragma region SETTERS/GETTERS
  
  void Problem::allgather(){
    
    Vec vout;
    VecScatter vecscat;
    
    VecScatterCreateToAll(p_vec[0],&vecscat,&vout);
    VecScatterBegin(vecscat,p_vec[0],vout,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(vecscat,p_vec[0],vout,INSERT_VALUES,SCATTER_FORWARD);
    
    VecScatterDestroy(&vecscat);
    VecCopy(vout,m_vecSeqTmp);
    VecDestroy(&vout);
    
  }
  
  
  double Problem::getTime(){
    return m_time;
  }
  
  double Problem::getTimeStep(){
    return m_dt;
  }
  
  void Problem::setOrderBDf(int order){
    m_orderbdf=order;
  }
  
  void Problem::setTimeStep(double dt){
    m_dt=dt;
  }
  
  void Problem::setNmax(int Nmax){
    m_Nmax=Nmax;
  }
  
  void Problem::setTmax(double Tmax){
    m_Tmax=Tmax;
  }
  
  void Problem::setProblemName(string pb_name){
    m_name = pb_name;
  }
  
#pragma endregion
  
  //#####################//
  //#####################//
  
#pragma region OHTERS
  void Problem::readInitComp(){
    
    m_comp1Init.clear();

    stringstream fileName;
    //fileName << "initVm_p" << m_pythonProc << ".txt"; 
    fileName << "./hetero/initVm_p" << m_pythonProc << ".txt"; 
    ifstream inputFile(fileName.str().c_str());
    if (!inputFile){
      cout << "Error. File " << fileName.str() << " does not exist." <<endl;
      exit(1);
    }
    
    while (true){
      double val;
      inputFile >> val;
      m_comp1Init.push_back(val);
      if (inputFile.eof()) break;
    }
    inputFile.close();    
      
  }
#pragma endregion
  
  //#####################//
  //#####################//
  
  
 
}
