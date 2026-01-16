#ifndef BETHE_BLOCH_H
#define BETHE_BLOCH_H

#include <TMath.h>

namespace dedx_constants
{
  // hadron masses
  constexpr double m_pi = 0.1396; // GeV
  constexpr double m_K = 0.4937; // GeV
  constexpr double m_p = 0.9382; // GeV
  constexpr double m_d = 1.876; // GeV

  // electron mass [eV]
  constexpr double m_e = 511e3;

  // TPC gas fractions
  constexpr double ar_frac = 0.75;
  constexpr double cf4_frac = 0.2;
  constexpr double isobutane_frac = 0.05;

  // Mean excitation [src: W. Blum, W. Riegler, L. Rolandi, "Particle Detection with Drift Chambers"]
  constexpr double ar_I = 188; // eV
  constexpr double cf4_I = 115; // eV
  constexpr double isobutane_I = 48.3; // eV

  // Mean excitation of mixture approximated using Bragg additivity rule
  constexpr double sphenix_I = ar_frac*ar_I + cf4_frac*cf4_I + isobutane_frac*isobutane_I;
}

// Bethe-Bloch fit function, vs. betagamma
// A = normalization constant, equal to (ADC conversion)*4pi*n*Z^2*e^4/(m_e*c^2*4pi*epsilon_0^2)
// B = A*(ln(2*m_e/I)-1) - (zero-suppression loss factor)
inline const double bethe_bloch_new(const double betagamma, const double A, const double B, const double C)
{
  const double beta = betagamma/sqrt(1.+betagamma*betagamma);

  return A/(beta*beta)*TMath::Log(betagamma) + A/(beta*beta)*B - A - C;
}

inline const double bethe_bloch_new_2D(const double betagamma, const double A, const double B)
{
  const double beta = betagamma/sqrt(1.+betagamma*betagamma);

  return A/(beta*beta)*(2.*TMath::Log(2.*dedx_constants::m_e/dedx_constants::sphenix_I * betagamma) - beta*beta) - B;
}

inline const double bethe_bloch_new_1D(const double betagamma, const double A)
{
  const double beta = betagamma/sqrt(1.+betagamma*betagamma);

  return A/(beta*beta)*(2.*TMath::Log(2.*dedx_constants::m_e/dedx_constants::sphenix_I * betagamma) - beta*beta);
}

// dE/dx for one gas species, up to normalization
inline const double bethe_bloch_species(const double betagamma, const double I)
{
  const double m_e = 511e3; // eV

  const double beta = betagamma/sqrt(1.+betagamma*betagamma);

  return 1./(beta*beta)*(TMath::Log(2.*m_e/I*betagamma*betagamma)-beta*beta);
}

// dE/dx for TPC gas mixture, up to normalization
inline const double bethe_bloch_total(const double betagamma)
{
  return dedx_constants::ar_frac * bethe_bloch_species(betagamma,dedx_constants::ar_I) +
         dedx_constants::cf4_frac * bethe_bloch_species(betagamma,dedx_constants::cf4_I) +
         dedx_constants::isobutane_frac * bethe_bloch_species(betagamma,dedx_constants::isobutane_I);
}

inline Double_t bethe_bloch_new_wrapper(Double_t* x, Double_t* par)
{
  Double_t betagamma = x[0];
  Double_t A = par[0];
  Double_t B = par[1];
  Double_t C = par[2];

  return bethe_bloch_new(betagamma,A,B,C);
}

inline Double_t bethe_bloch_new_2D_wrapper(Double_t* x, Double_t* par)
{
  Double_t betagamma = x[0];
  Double_t A = par[0];
  Double_t B = par[1];

  return bethe_bloch_new_2D(betagamma,A,B);
}

inline Double_t bethe_bloch_new_1D_wrapper(Double_t* x, Double_t* par)
{
  Double_t betagamma = x[0];
  Double_t A = par[0];

  return bethe_bloch_new_1D(betagamma,A);
}

// wrapper function for TF1 constructor, for fitting
inline Double_t bethe_bloch_wrapper(Double_t* ln_bg, Double_t* par)
{
  Double_t betagamma = exp(ln_bg[0]);

  Double_t norm = par[0];

  return norm * bethe_bloch_total(betagamma);
}

inline Double_t bethe_bloch_vs_p_wrapper(Double_t* x, Double_t* par)
{
  Double_t p = x[0];
  Double_t norm = par[0];
  Double_t m = par[1];

  return norm * bethe_bloch_total(fabs(p)/m);
}

inline Double_t bethe_bloch_vs_logp_wrapper(Double_t* x, Double_t* par)
{
  Double_t p = pow(10.,x[0]);
  Double_t norm = par[0];
  Double_t m = par[1];

  return norm * bethe_bloch_total(fabs(p)/m);
}

inline Double_t bethe_bloch_vs_p_wrapper_ZS(Double_t* x, Double_t* par)
{
  Double_t p = x[0];
  Double_t norm = par[0];
  Double_t m = par[1];
  Double_t ZS_loss = par[2];

  return norm * bethe_bloch_total(fabs(p)/m) - ZS_loss;
}

inline Double_t bethe_bloch_vs_p_wrapper_new(Double_t* x, Double_t* par)
{
  Double_t p = x[0];
  Double_t A = par[0];
  Double_t B = par[1];
  Double_t C = par[2];
  Double_t m = par[3];

  return bethe_bloch_new(fabs(p)/m,A,B,C);
}

inline Double_t bethe_bloch_vs_p_wrapper_new_2D(Double_t* x, Double_t* par)
{
  Double_t p = x[0];
  Double_t A = par[0];
  Double_t B = par[1];
  Double_t m = par[2];

  return bethe_bloch_new_2D(fabs(p)/m,A,B);
}

inline Double_t bethe_bloch_vs_p_wrapper_new_1D(Double_t* x, Double_t* par)
{
  Double_t p = x[0];
  Double_t A = par[0];
  Double_t m = par[1];

  return bethe_bloch_new_1D(fabs(p)/m,A);
}

inline Double_t bethe_bloch_vs_logp_wrapper_ZS(Double_t* x, Double_t* par)
{
  Double_t p = pow(10.,x[0]);
  Double_t norm = par[0];
  Double_t m = par[1];
  Double_t ZS_loss = par[2];

  return norm * bethe_bloch_total(fabs(p)/m) - ZS_loss;
}

inline Double_t bethe_bloch_vs_logp_wrapper_new(Double_t* x, Double_t* par)
{
  Double_t p = pow(10.,x[0]);
  Double_t A = par[0];
  Double_t B = par[1];
  Double_t C = par[2];
  Double_t m = par[3];

  return bethe_bloch_new(fabs(p)/m,A,B,C);
}

inline Double_t bethe_bloch_vs_logp_wrapper_new_1D(Double_t* x, Double_t* par)
{
  Double_t p = pow(10.,x[0]);
  Double_t A = par[0];
  Double_t m = par[1];

  return bethe_bloch_new_1D(fabs(p)/m,A);
}

// ratio of dE/dx between two particle species at the same momentum
// (useful for dE/dx peak fits)
inline const double dedx_ratio(const double p, const double m1, const double m2)
{
  const double betagamma1 = fabs(p)/m1;
  const double betagamma2 = fabs(p)/m2;

  return bethe_bloch_total(betagamma1)/bethe_bloch_total(betagamma2);
}

#endif // BETHE_BLOCH_H
