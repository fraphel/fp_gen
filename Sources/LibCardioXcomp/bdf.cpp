#include "stdafx.h"
#include "bdf.hpp"

namespace cardioxcomp {


  //#############//
  // CONSTRUCTOR //
  //#############//

  Bdf::Bdf():
    m_order(0),
    m_numberOfComp(1),
    m_size(0),
    m_build_order_2(false),
    m_build_order_3(false),
    m_build_RHS(false) {
  
  }

  //#############//
  //#############//


  //############//
  // DESTRUCTOR //
  //############//

  Bdf::~Bdf() {
    m_alpha.clear();
    m_beta.clear();
    for (int i=0 ; i< m_numberOfComp; i++) {
      if (m_build_order_2)
        VecDestroy(&m_sol_n_1[i]);
      if (m_build_order_3)
        VecDestroy(&m_sol_n_2[i]);
      if (m_build_RHS) {
        VecDestroy(&m_sol_n[i]);
        VecDestroy(&m_rhs[i]);
      }
    }
  }

  //############//
  //############//


  //###########//
  // BDF ORDER //
  //###########//

  void Bdf::defineOrder(int n, int nComp) {
    m_order = n;
    m_numberOfComp=nComp;
    m_alpha.resize(n+1);
    m_beta.resize(n);
    m_sol_n.resize(nComp);
    m_sol_n_1.resize(nComp);
    m_sol_n_2.resize(nComp);
    m_rhs.resize(nComp);
  
    switch (n) {
    case 1:
      m_alpha[0] = 1.0; // Backward Euler
      m_alpha[1] = 1.0;
      m_beta[0] = 1.0;// u^{n+1} \approx u^n
      break;
    case 2:
      m_alpha[0] = 3.0 / 2.0;
      m_alpha[1] = 2.0;
      m_alpha[2] = -1.0 / 2.0;
      m_beta[0] = 2.0;
      m_beta[1] = -1.0;
      break;
    case 3:
      m_alpha[0] = 11.0 / 6.0;
      m_alpha[1] = 3.0;
      m_alpha[2] = -3.0 / 2.0;
      m_alpha[3] = 1.0 / 3.0;
      m_beta[0] = 3.0;
      m_beta[1] = -3.0;
      m_beta[2] = 1.0;
      break;
    }
  }

  //###########//
  //###########//


  //################//
  // INITIALIZE BDF //
  //################//

  void Bdf::initialize(Vec& sol_0) {
    switch (m_order) {
    case 1:
      VecGetSize(sol_0,&m_size);

      VecDuplicate(sol_0,&m_sol_n[0]);
      VecCopy(sol_0,m_sol_n[0]); 

      VecDuplicate(sol_0,&m_rhs[0]);

      m_build_RHS = true;
      break;
    case 2:
      VecGetSize(sol_0,&m_size);

      VecDuplicate(sol_0,&m_sol_n[0]);
      VecCopy(sol_0,m_sol_n[0]); 

      VecDuplicate(sol_0,&m_sol_n_1[0]);
      VecCopy(sol_0,m_sol_n_1[0]); 
      m_build_order_2 = true;

      VecDuplicate(sol_0,&m_rhs[0]);
      m_build_RHS = true;
      break;
    case 3:
      VecGetSize(sol_0,&m_size);

      VecDuplicate(sol_0,&m_sol_n[0]);
      VecCopy(sol_0,m_sol_n[0]); 

      VecDuplicate(sol_0,&m_sol_n_1[0]);
      VecCopy(sol_0,m_sol_n_1[0]); 
      m_build_order_2 = true;

      VecDuplicate(sol_0,&m_sol_n_2[0]);
      VecCopy(sol_0,m_sol_n_2[0]); 
      m_build_order_3 = true;

      VecDuplicate(sol_0,&m_rhs[0]);
      m_build_RHS = true;

      break;
    }
  }

  //################//
  //################//


  //##################//
  // UPDATE NEW VALUE //
  //##################//

  void Bdf::update(Vec& new_sol_n) {

    if(m_build_order_3) {
      VecCopy(m_sol_n_1[0],m_sol_n_2[0]);
      VecCopy(m_sol_n[0],m_sol_n_1[0]);
      VecCopy(new_sol_n,m_sol_n[0]);

    } else if(m_build_order_2) {

      VecCopy(m_sol_n[0],m_sol_n_1[0]);
      VecCopy(new_sol_n,m_sol_n[0]);

    } else if(m_build_RHS) {

      VecCopy(new_sol_n,m_sol_n[0]);

    }
  }

  //##################//
  //##################//


  //#############//
  // COMPUTE RHS //
  //#############//

  void Bdf::computeRHSTime(double dt, Vec& RHSTime) {

    VecZeroEntries(RHSTime);

    switch (m_order) {
    case 1:
      VecAXPY(RHSTime,m_alpha[1]/dt,m_sol_n[0]); // RHSTime = RHSTime + (m_alpha/dt)*m_sol_n[0] 
      break;
    case 2:
      VecAXPY(RHSTime,m_alpha[1]/dt,m_sol_n[0]);
      VecAXPY(RHSTime,m_alpha[2]/dt,m_sol_n_1[0]);
      break;
    case 3:
      VecAXPY(RHSTime,m_alpha[1]/dt,m_sol_n[0]);
      VecAXPY(RHSTime,m_alpha[2]/dt,m_sol_n_1[0]);
      VecAXPY(RHSTime,m_alpha[3]/dt,m_sol_n_2[0]);
      break;
    }
  }


  void Bdf::computeRHSTime(double dt) {
    computeRHSTime(dt,m_rhs[0]);
  }

  //#############//
  //#############//


  //###################//
  // EXTRAPOLATE VALUE //
  //###################//

  void Bdf::extrapolate(Vec& extrap) {
    VecZeroEntries(extrap);

    switch (m_order) {
    case 1:
      VecAXPY(extrap,m_beta[0],m_sol_n[0]);  // extrap=extrap+m_sol_n[0]*m_beta[0] (bdf1:beta[0]=1)
      break;
    case 2:
      VecAXPY(extrap,m_beta[0],m_sol_n[0]);
      VecAXPY(extrap,m_beta[1],m_sol_n_1[0]);
      break;
    case 3:
      VecAXPY(extrap,m_beta[0],m_sol_n[0]);
      VecAXPY(extrap,m_beta[1],m_sol_n_1[0]);
      VecAXPY(extrap,m_beta[2],m_sol_n_2[0]);
      break;
    }
  }

  //###################//
  //###################//


}
