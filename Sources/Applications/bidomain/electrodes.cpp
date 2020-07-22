#include "stdafx.h"

#include "electrodes.hpp"


namespace cardioxcomp{

  //#############//
  // CONSTRUCTOR //
  //#############//

  Electrodes::Electrodes(){

  }

  //#############//
  //#############//


  //############//
  // DESTRUCTOR //
  //############//

  Electrodes::~Electrodes(){

    VecDestroy(&m_Iel_Extrap);
    VecDestroy(&m_Ue_Extrap);
    VecDestroy(&m_Iel);
    VecDestroy(&m_Ue);
    VecDestroy(&dU);

  }

  //############//
  //############//


  //################//
  // INITIALIZATION //
  //################//

  void Electrodes::initializeElectrodes(int nE, double dt, const vector<double>& C, const vector<double>& Ri, const vector<double>& Rel, const vector<double>& diameters){
    
    m_nbElectrodes = nE;
    m_dt = dt;
    
    m_C = C;
    m_Rel = Rel;
    m_Ri = Ri;
    
    
    m_diameters = diameters;
    computeAreas();
    computeTau();
    computeRtilde();
    
    
    m_dU.resize(m_nbElectrodes);
    m_Umean.resize(m_nbElectrodes);
    
    for(int ie=0 ; ie<m_nbElectrodes ; ie++){
      m_dU[ie] = 0.;
      m_Umean[ie] = 0.;
    }
    
    // CREATION DES VECTEURS
    //_________________________________________________________________
    
    //Pour Ue
    VecCreate(PETSC_COMM_SELF,&m_Ue);
    VecSetSizes(m_Ue,PETSC_DECIDE,m_nbElectrodes);
    
    VecSetFromOptions(m_Ue);
    VecSet(m_Ue,0.);
    
    VecAssemblyBegin(m_Ue);
    VecAssemblyEnd(m_Ue);
    
    
    VecCreate(PETSC_COMM_SELF,&m_Ue_Extrap);
    VecSetSizes(m_Ue_Extrap,PETSC_DECIDE,m_nbElectrodes);
    
    VecSetFromOptions(m_Ue_Extrap);   
    VecSet(m_Ue_Extrap,0.);
    
    VecAssemblyBegin(m_Ue_Extrap);
    VecAssemblyEnd(m_Ue_Extrap);
    
    //Pour Iel
    VecCreate(PETSC_COMM_SELF,&m_Iel);
    VecSetSizes(m_Iel,PETSC_DECIDE,m_nbElectrodes);
    
    VecSetFromOptions(m_Iel);
    VecSet(m_Iel,0.);
    
    VecAssemblyBegin(m_Iel);
    VecAssemblyEnd(m_Iel);
    
    
    VecCreate(PETSC_COMM_SELF,&m_Iel_Extrap);
    VecSetSizes(m_Iel_Extrap,PETSC_DECIDE,m_nbElectrodes);
    
    VecSetFromOptions(m_Iel_Extrap);   
    VecSet(m_Iel_Extrap,0.);
    
    VecAssemblyBegin(m_Iel_Extrap);
    VecAssemblyEnd(m_Iel_Extrap);
    
    //For dU
    VecCreate(PETSC_COMM_SELF,&dU);
    VecSetSizes(dU,PETSC_DECIDE,m_nbElectrodes);
    
    VecSetFromOptions(dU);
    VecSet(dU,0.);
    
    VecAssemblyBegin(dU);
    VecAssemblyEnd(dU);
    
    //_________________________________________________________________
    //_________________________________________________________________
    
    
    m_bdfIel.defineOrder(3,1);
    m_bdfUe.defineOrder(3,1);
    
    //For Ue
    m_bdfUe.initialize(m_Ue);
    m_bdfUe.extrapolate(m_Ue_Extrap);
    
    //For Iel    
    m_bdfIel.initialize(m_Iel);
    m_bdfIel.extrapolate(m_Iel_Extrap);
    
    
    //---------------------------------------------------------------


    
  }

  //################//
  //################//


  //###############//
  // COMPUTE AREAS //
  //###############//

  void Electrodes::computeAreas(){

    m_areas.resize(m_nbElectrodes);

    for(int ek=0 ; ek<m_nbElectrodes ; ek++){
      const double rk = 0.5*m_diameters[ek];
      m_areas[ek] = rk*rk*4.*atan(1.);
    }
  }

  //###############//
  //###############//


  //#########################//
  // COMPUTE TAU = C(Ri+Rel) //
  //#########################//

  void Electrodes::computeTau(){

    m_tauk.resize(m_nbElectrodes);

    for(int ek=0 ; ek<m_nbElectrodes ; ek++){
      const double tauk = (m_Rel[ek]+m_Ri[ek])*m_C[ek];
      m_tauk[ek] = tauk;
    }
  }

  //#########################//
  //#########################//


  //#########################//
  // COMPUTE RTILDE = Ri+Rel //
  //#########################//

  void Electrodes::computeRtilde(){

    m_rtilde.resize(m_nbElectrodes);
    for(int ek=0 ; ek<m_nbElectrodes ; ek++){
      const double rtilde = m_Rel[ek]+m_Ri[ek];
      m_rtilde[ek] = rtilde;
    }

  }

  //#########################//
  //#########################//


  //##########################//
  // FORWARD: ADVANCE IN TIME //
  //##########################//

  void Electrodes::forwardElectrodes(){
    
    // Updating and extrapolation
    m_bdfIel.update(m_Iel);
    m_bdfIel.extrapolate(m_Iel_Extrap);
    
    
    m_bdfUe.update(m_Ue); // time derivative of Ue
    m_bdfUe.extrapolate(m_Ue_Extrap);
    //__________________________
    
    
    // dU calculation
    computedU();
    //__________________________
    
    
    // Iel calculation
    computeIel();
    //__________________________
    
    
    //Calcul de Umes
    computeUmes();
    //__________________________
    
  }
  
  //##########################//
  //##########################//
  
  
  //#############//
  // Compute Iel //
  //#############//

  void Electrodes::computeIel(){

    platformInt locsize;
    platformInt minnum,maxnum;

    VecGetLocalSize(m_bdfIel.rhs(),&locsize);
    VecGetOwnershipRange(m_bdfIel.rhs(),&minnum,&maxnum);

    //~~~~~~~~~

    m_bdfIel.computeRHSTime(m_dt);

    VecAXPY(m_bdfIel.rhs(),1,dU);

    for(int ie=0 ; ie<locsize ; ie++){
    
      double& coeffDeriv = m_bdfIel.coeffDeriv0();

      platformInt pos = minnum+ie;
      double value_RHS;
      VecGetValues(m_bdfIel.rhs(),1,&pos,&value_RHS);
      
      const double coeff = (coeffDeriv/m_dt)+1./m_tauk[ie];

      //const double Cverre = 10e-9;
      //const double coeff = (coeffDeriv/m_dt)+(m_C[ie]+Cverre)/(m_C[ie]*Cverre*(m_Rel[ie]+m_Ri[ie]));

      const double Iel_np1 = value_RHS/coeff;
      VecSetValues(m_Iel,1,&pos,&Iel_np1, INSERT_VALUES);

    }
    
    VecAssemblyBegin(m_Iel);
    VecAssemblyEnd(m_Iel);
    
  }

  //#############//
  //#############//


  //##########################//
  // COMPUTE dU : FOR RHS IEL //
  //##########################//

  void Electrodes::computedU(){

    m_bdfUe.computeRHSTime(m_dt);
    
    platformInt locsize;
    platformInt minnum,maxnum;

    VecGetLocalSize(m_bdfUe.rhs(),&locsize);
    VecGetOwnershipRange(m_bdfUe.rhs(),&minnum,&maxnum);
    
    for(int ek=0 ; ek<locsize ; ek++){

      double& coeffDeriv = m_bdfUe.coeffDeriv0();
      platformInt pos = minnum+ek;
      double value_RHS;
      VecGetValues(m_bdfUe.rhs(),1,&pos,&value_RHS);

      double Um = m_Umean[ek];
      m_dU[ek] = Um*(coeffDeriv/m_dt)-value_RHS;
      VecSetValue(dU, pos, m_dU[ek]/m_rtilde[ek], INSERT_VALUES);
    
    }

    VecAssemblyBegin(dU);
    VecAssemblyEnd(dU);

  }

  //##########################//
  //##########################//


  //##############//
  // COMPUTE UMES //
  //##############//

  void Electrodes::computeUmes(){

    m_Umes.clear();
    
    platformInt locsize;
    platformInt minnum,maxnum;

    VecGetLocalSize(m_Iel_Extrap,&locsize);
    VecGetOwnershipRange(m_Iel_Extrap,&minnum,&maxnum);

    for(int ie=0 ; ie<locsize ; ie++){

      platformInt pos = minnum+ie;
      double Iel;
      VecGetValues(m_Iel_Extrap,1,&pos,&Iel);

      m_Umes.push_back(m_Ri[ie]*Iel);

      double val = m_Umean[ie];
      VecSetValues(m_Ue,1,&pos,&val, INSERT_VALUES);

    }

    VecAssemblyBegin(m_Ue);
    VecAssemblyEnd(m_Ue);

  }

  //##############//
  //##############//


  //#####################//
  // COMPUTE SOURCE TERM //
  //#####################//

  void Electrodes::computeSourceIel(vector<double>& source){

    source.clear();
    
    platformInt locsize;
    platformInt minnum,maxnum;

    VecGetLocalSize(m_Iel_Extrap,&locsize);
    VecGetOwnershipRange(m_Iel_Extrap,&minnum,&maxnum);

    for(int ie=0 ; ie<locsize ; ie++){
      platformInt pos = minnum+ie;
      double Iel;
      VecGetValues(m_Iel_Extrap,1,&pos,&Iel);

      source.push_back(Iel/m_areas[ie]);
    }

  }

  //#####################//
  //#####################//


  //#####################//
  // GETTERS AND SETTERS //
  //#####################//

  void Electrodes::getUmes(vector<double>& Umes){

    Umes = m_Umes;

  }

  void Electrodes::set_Umean(const vector<double>& umean){

    m_Umean = umean;

  }

  //#####################//
  //#####################//




}
