#ifndef _BDF_HPP
#define _BDF_HPP

//or #pragma once and no #endif at the end!

#include <petscvec.h>
#include <petscsys.h>

#include<vector>
#include<iostream>
#include "platform.h"

using namespace std;

const int BDF_MAX_ORDER = 3;

namespace cardioxcomp {

class Bdf {
public:
  Bdf();
  ~Bdf();

  void defineOrder(int n, int nComp=1);
  void initialize(Vec& sol_0);


  void update(Vec& sol_n);  
  void computeRHSTime(double dt, Vec& RHSTime);
  void computeRHSTime(double dt);
  void extrapolate(Vec& extrap);


   
  //##################//
  // Access functions //
  //##################//

 
  //return \alpha_0.
  inline const double & coeffDeriv0() const {
    return m_alpha[0];
  }
  inline double & coeffDeriv0() {
    return m_alpha[0];
  }
  
  // return alpha[i]
  inline const double & alpha(int i) const {
    return m_alpha[i];
  }
  inline double & alpha(int i) {
    return m_alpha[i];
  }
  
  // return beta[i]
  inline const double & beta(int i) const {
    return m_beta[i];
  }
  inline double & beta(int i) {
    return m_beta[i];
  }

  inline const Vec & sol_n() const {
    return m_sol_n[0];
  }
  inline Vec & sol_n() {
    return m_sol_n[0];
  }
  
  inline const Vec & sol_n_1() const {
    return m_sol_n_1[0];
  }
  inline Vec & sol_n_1() {
    return m_sol_n_1[0];
  }
  
  inline const Vec & sol_n_2() const {
    return m_sol_n_2[0];
  }
  inline Vec & sol_n_2() {
    return m_sol_n_2[0];
  }
  
  inline const Vec & rhs() const {
    return m_rhs[0];
  }
  inline Vec & rhs() {
    return m_rhs[0];
  }
  
  inline const int & order() const {
    return m_order;
  }
  inline int & order() {
    return m_order;
  }
  
  inline const int & numComp() const {
    return m_numberOfComp;
  }
  inline int & numComp() {
    return m_numberOfComp;
  }

  //##################//
  //##################//

private:
  //! Order of the BDF derivative/extrapolation: the time-derivative
  //! coefficients vector has size n+1, the extrapolation vector has size n
  int m_order;
  //! Number of components
  int m_numberOfComp;
  //! Size du vecteur solution
  platformInt m_size;
  //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
  std::vector<double> m_alpha;
  //! Coefficients \f$ \beta_i \f$ of the extrapolation
  std::vector<double> m_beta;
  
  std::vector<Vec> m_sol_n;
  std::vector<Vec> m_sol_n_1;
  std::vector<Vec> m_sol_n_2;
  std::vector<Vec> m_rhs;

  bool m_build_order_2;
  bool m_build_order_3;
  bool m_build_RHS;
};


}


#endif
