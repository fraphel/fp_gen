#ifndef __CardioXcomp__problem__
#define __CardioXcomp__problem__

#ifdef _WIN32
#pragma warning(disable : 4996)
#endif

#include <iostream>
#include<iomanip>

#include <vector>

#include <petscsys.h>
#include <petscksp.h>

#include "mesh.hpp"
#include "dof.hpp"
#include "matElem.hpp"
#include "vecElem.hpp"
#include "boundaryCondition.hpp"
#include <sstream>
#include "bdf.hpp"

#include "output.hpp"

#include "platform.h"

#include "luafile.hpp"

#include <mpi.h>

using namespace std;

namespace cardioxcomp {
  
  class Problem{
  protected:
    Mesh* p_mesh;  // class Mesh
    Mesh* p_localMesh;
    Dof* p_dof;
    Mat* p_mat;
    Vec* p_vec;
    Vec m_vecSeqTmp;
    KSP m_ksp;
    PC m_pc;
    int m_numComp;
    int m_numMat;
    int m_numVec;
    int m_numProc;
    int m_rankProc;
    MatElem* p_matElem;
    VecElem* p_vecElem;
    std::vector<PetscInt> m_loc2GlobElem; // useful ?
    //! mapping between local to global ordering of element.
    ISLocalToGlobalMapping m_mappingElem;
    bool m_initMappingElem;

    //Post-treatment variables
    std::vector<string> m_physicalvariables;
    int m_freq_writeSolution;

    //Time variables
    Bdf m_bdfEdp;
    int m_orderbdf;
    double m_dt;
    Vec m_extrapolate;

    Vec m_source;
    double m_time;
    int m_numIt;
    double m_Tmax;
    int m_Nmax;
    string m_name;
    
    std::vector<double> m_physicalInitialValues;
    std::vector<double> m_heteroInitial;
    std::vector<double> m_comp1Init;
    
    double m_coefftime,m_coeffsource;
    
    // Solver
    //____________
    
    double m_edpReltol,m_edpAbstol;
    int m_edpMaxiter;
    string m_solver,m_preconditioner,m_preconditionerOption;
    MatStructure m_setPreconditionerOption;
    
    //____________
    //____________
    
    // Parallel
    //__________
    
    platformInt m_localsize;
    platformInt m_minOwnership;
    string m_pythonProc;
    
    //__________
    //__________
    
    // Theta
    double m_theta;
    bool m_isHetero;
    bool m_cellHetero;
    
    
  public:
    Problem();
    Problem(std::string meshfile,int numComp,int numMat,int numVec,int numProc,int rankProc,bool isHetero, string pythonProc, bool cellHetero);

#ifdef _WIN32
    __declspec(dllexport) 
#endif
    ~Problem();

    virtual void computeMatVec() = 0;
    virtual void dirichletBC() = 0;
    
    virtual Vec computeODEpart(double,Vec&,bool) = 0;
    virtual void initializeProblem(const Vec&) = 0;
    virtual void writeElectrodesSolution() = 0;
    virtual void writeElectrodesSolution_headers(){};
    virtual void writeUemesAtTheEnd(){};

    virtual void setLuaConfig(const LuaFile&) = 0;
    virtual void setLuaVar(const LuaVariables&) = 0;
    
    void removeMean(int comp=1);
    virtual double computeCaiMean() = 0;
    
    std::vector< std::vector<double> > SolutionElectrodes;
    std::vector< std::vector<double> > SolutionVm;
    
    void buildSolver();
    void finalizeMatVec();
    void solve();
    void postprocScalar(const std::string& filename,int ivec,int icomp);

    void forward();
    void forward_split();

    void setOrderBDf(int order);
    void setTimeStep(double dt);
    void setNmax(int Nmax);
    void setTmax(double Tmax);
    void setProblemName(string pb_name);
    
    void allgather();
        
    double getTime();
    int m_simID;
    void setTime(double newTime){m_time=newTime;};
    void setNumIt(int newIt){m_numIt = newIt;};
    double getTimeStep();

    void readInitComp();
    
  };
  
}

#endif
