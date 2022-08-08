
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

#ifndef EVTVUBAC_HH
#define EVTVUBAC_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

#include <vector>

class EvtParticle;

// Analytic Coupling Model (based on hep-ph/0608047 by Aglietti, Ferrera and Ricciardi)

class EvtVubAC : public EvtDecayIncoherent {
  public:
    std::string getName() override;

    EvtDecayBase* clone() override;

    void initProbMax() override;

    void init() override;

    void decay( EvtParticle* Bmeson ) override;

  private:
    // Input parameters
    double mB;

    double alphaSmZ;
    double alphaSmB;
    double c;
    double q;
    double k;

    double CF;
    double CA;

    double beta0;

    std::vector<double> gvars;

    double rate( double u, double w, double xb );
    double wreg( double w );
    double alphaS( double Q );
    double PolyLog( double v, double z );
    double ureg( double u );
    double ularge( double u );
    double Coeff( double u, double w, double xb );
    double Coeff1( double w, double xb );
    double Coeff0( double w, double xb );
    double Sigma( double x1, double x2 );
    double max( double ub, double lb );
    double d1( double u, double w, double xb );
    double d( double u, double w, double xb );
    double f( double w );
    double Lambda2( double x, double alphaSmZ );
    int Bisect( double x1, double x2, double precision, double& root,
                const double alphaSmZ );
    double FindRoot( const double alphaSmZ );
};

#endif
