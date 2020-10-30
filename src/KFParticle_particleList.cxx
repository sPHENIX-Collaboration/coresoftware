/*
 * This file is part of KFParticle package
 * Copyright ( C ) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * ( at your option ) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/*****************/
/* Cameron Dean  */
/*   LANL 2020   */
/* cdean@bnl.gov */ 
/*****************/

/** Mass Hypothesis Codes
 ** These are the codes used in the PDG to identify different particles
 ** that have also been coded into KFParticle
 ** PDG Code | Mass Index | Particle
 **         11 |  0 | Electron
 **         13 |  1 | Muon
 **         19 |  1 | Muon
 **        211 |  2 | Pion ( Charged )
 **        321 |  3 | Kaon ( Charged )
 **       2212 |  4 | Proton
 ** 1000010020 |  5 | Deuteron
 ** 1000010030 |  6 | Triton
 ** 1000020030 |  7 | Helium-3
 ** 1000020040 |  8 | Helium-4
 **       3112 |  9 | Sigma ( - )
 **       3222 | 10 | Sigma ( + )
 **       3312 | 11 | Xi
 **       3334 | 12 | Omega
 ** Use PDG codes in your analysis
 ** All other codes return a pion
 **/

#include "KFParticle_particleList.h"

//KFParticle stuff
#include "KFParticleDatabase.h"

/// KFParticle constructor
KFParticle_particleList::KFParticle_particleList(){}

/// KFParticle destructor
KFParticle_particleList::~KFParticle_particleList(){}

std::map<std::string, float> KFParticle_particleList::getParticleList()
{
  //There is a scope issue here and kfpDatabase must be declared in function unlike in KFPTools
  KFParticleDatabase kfpDatabase;
  std::map<std::string, float> particleMasses; 

  //Leptons
  particleMasses["electron"] = kfpDatabase.GetMass(11);
  particleMasses["muon"]     = kfpDatabase.GetMass(13);
  particleMasses["tau"]      = 1.77686;

  //Gauge bosons and Higgs
  particleMasses["W+"]    = 80.379;
  particleMasses["W-"]    = 80.379;
  particleMasses["Z"]     = 91.1876;
  particleMasses["Higgs"] = 125.10;

  //Light, unflavoured mesons
  particleMasses["pion"]      = kfpDatabase.GetMass(211);
  particleMasses["pi+"]       = kfpDatabase.GetMass(211);
  particleMasses["pi-"]       = kfpDatabase.GetMass(211);
  particleMasses["pi0"]       = kfpDatabase.GetPi0Mass();
  particleMasses["eta"]       = 0.547862;
  particleMasses["f0(500)"]   = 0.5;
  particleMasses["rho"]       = 0.77526;
  particleMasses["rho(770)"]  = 0.77526;
  particleMasses["f0(980)"]   = 0.990;
  particleMasses["phi"]       = 1.019461;
  particleMasses["phi(1020)"] = 1.019461;

  //Strange mesons   
  particleMasses["kaon"]    = kfpDatabase.GetMass(321);
  particleMasses["K+"]      = kfpDatabase.GetMass(321);
  particleMasses["K-"]      = kfpDatabase.GetMass(321);
  particleMasses["K0"]      = 0.497611;
  particleMasses["KS0"]     = 0.497611;
  particleMasses["KL0"]     = 0.497611;
  particleMasses["K*(892)"] = 0.89166;

  //Light baryons
  particleMasses["proton"]  = kfpDatabase.GetMass(2212);
  particleMasses["neutron"] = 0.93957;
  particleMasses["Lambda"]  = 1.11568;
  particleMasses["Sigma+"]  = kfpDatabase.GetMass(3222);
  particleMasses["Sigma0"]  = 1.192642;
  particleMasses["Sigma-"]  = kfpDatabase.GetMass(3112);
  particleMasses["Xi0"]     = 1.31486;
  particleMasses["Xi+"]     = 1.32171;
  particleMasses["Xi-"]     = 1.32171;

  //Charm-hadrons
  particleMasses["D0"]       = 1.86483;
  particleMasses["D0bar"]    = 1.86483;
  particleMasses["D+"]       = 1.86965;
  particleMasses["D-"]       = 1.86965;
  particleMasses["Ds+"]      = 1.96834;
  particleMasses["Ds-"]      = 1.96834;
  particleMasses["D*0"]      = 2.00685;
  particleMasses["D*+"]      = 2.01026;
  particleMasses["D*-"]      = 2.01026;
  particleMasses["Ds*+"]     = 2.11220;
  particleMasses["Ds*-"]     = 2.11220;
  particleMasses["Lc+"]      = 2.28646;
  particleMasses["Lambdac+"] = 2.28646;
  particleMasses["Xic0"]     = 2.47090;
  particleMasses["Xic+"]     = 2.46794;
  particleMasses["Xic-"]     = 2.46794;
  particleMasses["Omegac"]   = 2.6952;
  particleMasses["Xicc++"]   = 3.6212;

  //B-hadrons
  particleMasses["B+"]       = 5.279; 
  particleMasses["B-"]       = 5.279;
  particleMasses["B0"]       = 5.279;
  particleMasses["Bs0"]      = 5.366;
  particleMasses["Bc+"]      = 6.2749;
  particleMasses["Bc-"]      = 6.2749;
  particleMasses["Bc"]       = 6.2749;
  particleMasses["Bc(2S)"]   = 6.8716;
  particleMasses["Lambdab0"] = 5.61960;
  particleMasses["Sigmab+"]  = 5.81056;
  particleMasses["Sigmab-"]  = 5.81056;
  particleMasses["Xib0"]     = 5.7919;
  particleMasses["Xib+"]     = 5.7970;
  particleMasses["Xib-"]     = 5.7970;
  particleMasses["Omegab+"]  = 6.0461;
  particleMasses["Omegab-"]  = 6.0461;
  
  //Quarkonia
  //c-cbar
  particleMasses["J/psi"]       = 3.09690;
  particleMasses["psi(2S)"]     = 3.68610;
  particleMasses["X(3872)"]     = 3.87169;
  particleMasses["chic1(3872)"] = 3.87169;
  //b-bbar
  particleMasses["Upsilon(1S)"] = 9.46030;
  particleMasses["Upsilon(2S)"] = 10.02326;
  particleMasses["Upsilon(3S)"] = 10.3552;
  particleMasses["Upsilon(4S)"] = 10.5794;
  particleMasses["Upsilon(5S)"] = 10.8852;

  return particleMasses;

}


float KFParticle_particleList::returnPDGMass( const int pdgIndex ) ///Return mother masses from KFParticleDatabase ( pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf )
{ 
  KFParticleDatabase kfpDatabase;

   float mass, width;
   kfpDatabase.GetMotherMass( pdgIndex, mass, width );
   return mass;
 }
