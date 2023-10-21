#pragma once
#include"physicalconst.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>
#include <functional>
double B_nu(double nu) //  
{
    return (2*h*pow(nu,3)/pow(c,2))*(1/(exp(h*nu/(k*T)) - 1));
}
double J_nu(double r, double tau, double nu) 
{
   
     return 0.25*B_nu(nu)*pow(R_star/r,2)*exp(-tau);
    
}
double sigma(double nu) 
{
    return 6.8*1e-18*pow(v_ion/nu,3);
}
double sigmaS(double nu)//Verner 1996 
{
    double E=4.1356671*1e-15*nu;
    double E0=1.808*1e+1;
    double y0=9.935e-1;
   double sigma_0=4.564*1e+4;
   double x=E/E0-y0;
double y1=2.486e-1;
   double y=sqrt(x*x+y1*y1);
   double yw=6.385e-1;
   double ya=1.0;
   double P=1.361e+1;
    return sigma_0*(pow(x-1,2)+pow(yw,2))*pow(y,(0.5*P-5.5))*pow(1+sqrt(y/ya),-P)*1e-18;
}

double PHI(double r, double tau, double nu) 
{
    return 4*M_PI*(J_nu(r, tau, nu)/(h*nu))*sigma(nu);
    
}
double PHIS(double r, double tau, double nu) 
{
    return 4*M_PI*(J_nu(r, tau, nu)/(h*nu))*sigmaS(nu);
    
}
double alpha_B() 
{
    
    return 2.59e-13;// t_e nebula 10000;
}

double absorbt_func(double n_h,double sig,double r)
{
return n_h*sigma(sig);

}
 double _rrec[2][30][30];

  /*! @brief rnew array from Verner's script. */
  double _rnew[4][30][30];

  /*! @brief fe array from Verner's script. */
  double _fe[3][13];



void VernerRecombinationRates() {
  std::ifstream file("/home/anastasia/Downloads/CMacIonize-master/data/verner_rec_data.txt");

  std::string line;
  // skip comment line
  getline(file, line);

  // read rrec values
  {
    // skip comment line
    getline(file, line);
    for (uint_fast8_t i = 0; i < 2; ++i) {
      for (uint_fast8_t j = 0; j < 30; ++j) {
        getline(file, line);
        std::stringstream lstream(line);
        for (uint_fast8_t k = 0; k < 30; ++k) {
          lstream >> _rrec[i][j][k];
        }
      }
    }
  }

  // read rnew values
  {
    // skip comment line
    getline(file, line);
    for (uint_fast8_t i = 0; i < 4; ++i) {
      for (uint_fast8_t j = 0; j < 30; ++j) {
        getline(file, line);
        std::stringstream lstream(line);
        for (uint_fast8_t k = 0; k < 30; ++k) {
          lstream >> _rnew[i][j][k];
        }
      }
    }
  }

  // read fe values
  {
    // skip comment line
    getline(file, line);
    for (uint_fast8_t i = 0; i < 3; ++i) {
      getline(file, line);
      std::stringstream lstream(line);
      for (uint_fast8_t j = 0; j < 13; ++j) {
        lstream >> _fe[i][j];
      }
    }
  }

  // invert _rnew[2] and _rnew[3] values, as they are only used in divisions
  // and multiplying is much more efficient than dividing
  for (uint_fast8_t i = 0; i < 30; ++i) {
    for (uint_fast8_t j = 0; j < 30; ++j) {
      if (_rnew[2][i][j] != 0.) {
        _rnew[2][i][j] = 1. / _rnew[2][i][j];
      }
      if (_rnew[3][i][j] != 0.) {
        _rnew[3][i][j] = 1. / _rnew[3][i][j];
      }
    }
  }
}









double get_recombination_rate_verner(
    const uint_fast8_t iz, const uint_fast8_t in, const double T) 
     {

  double r = 0.;



  if (in <= 3 || in == 11 || (iz > 5 && iz < 9) || iz == 10 ||
      (iz == 26 && in > 11)) {
    const double tt = std::sqrt(T * _rnew[2][iz - 1][in - 1]);
    r = _rnew[0][iz - 1][in - 1] /
        (tt * std::pow(tt + 1., 1. - _rnew[1][iz - 1][in - 1]) *
         std::pow(1. + std::sqrt(T * _rnew[3][iz - 1][in - 1]),
                  1. + _rnew[1][iz - 1][in - 1]));
  } else {
    const double tt = T * 1.e-4;
    if (iz == 26 && in <= 13) {
      r = _fe[0][in - 1] *
          std::pow(tt, -_fe[1][in - 1] - _fe[2][in - 1] * std::log10(tt));
    } else {
      r = _rrec[0][iz - 1][in - 1] * std::pow(tt, -_rrec[1][iz - 1][in - 1]);
    }
  }

  return r;
}
double get_recombination_rate(
     const double temperature) 
      {
        const double T_in_eV = temperature / 1.16045221e4;
    double ratep1 = get_recombination_rate_verner(16, 15, temperature) +
           1.37e-9 * std::exp(-14.95 / T_in_eV) * std::pow(T_in_eV, -1.5);

    // Mazzotta et al. (1998), equation 7, validity range not specified
    
    double T_in_eV_inv = 1. / T_in_eV;
    double ratep2 = get_recombination_rate_verner(16, 14, temperature) +
           (8.0729e-9 * std::exp(-17.56 * T_in_eV_inv) +
            1.1012e-10 * std::exp(-7.07 * T_in_eV_inv)) *
               std::pow(T_in_eV, -1.5);

    // Abdel-Naby et al. (2012), equation 3, validity range [90 K; 9x10^7 K].
   
  
    // Abdel-Naby et al. (2012), equation 3, validity range [90 K; 9x10^7 K].
    double T_inv = 1. / temperature;
    
    double ratep3 = get_recombination_rate_verner(16, 13, temperature) +
           (5.817e-7 * std::exp(-362.8 * T_inv) +
            1.391e-6 * std::exp(-1058. * T_inv) +
            1.123e-5 * std::exp(-7160. * T_inv) +
            1.521e-4 * std::exp(-3.26e4 * T_inv) +
            1.875e-3 * std::exp(-1.235e5 * T_inv) +
            2.097e-2 * std::exp(-2.07e5 * T_inv)) *
               std::pow(temperature, -1.5);
  
  return std::max(0., ratep1+ratep2+ratep3);
      }
double alpha_B_Sulphur() 
{
    double delta = get_recombination_rate(10000);
    return delta;
}