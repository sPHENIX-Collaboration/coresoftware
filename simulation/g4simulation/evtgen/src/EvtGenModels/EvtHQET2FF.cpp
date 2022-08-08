
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

#include "EvtGenModels/EvtHQET2FF.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <math.h>
#include <stdlib.h>
#include <string>

EvtHQET2FF::EvtHQET2FF( double hqetrho2, double hqetha1_1, double hqetr1_1,
                        double hqetr2_1 )
{
    rho2 = hqetrho2;
    r1_1 = hqetr1_1;
    r2_1 = hqetr2_1;
    ha1_1 = hqetha1_1;
    extended = false;
}

EvtHQET2FF::EvtHQET2FF( double hqetrho2, double hqetha1_1, double hqetr1_1,
                        double hqetr2_1, double hqetr0_1 )
{
    rho2 = hqetrho2;
    r0_1 = hqetr0_1;
    r1_1 = hqetr1_1;
    r2_1 = hqetr2_1;
    ha1_1 = hqetha1_1;
    extended = true;
}

EvtHQET2FF::EvtHQET2FF( double hqetrho2, double hqetv1_1 )
{
    rho2 = hqetrho2;
    v1_1 = hqetv1_1;
    extended = false;
}

EvtHQET2FF::EvtHQET2FF( double hqetrho2, double hqetv1_1, double indelta )
{
    rho2 = hqetrho2;
    v1_1 = hqetv1_1;
    delta = indelta;
    extended = true;
}

void EvtHQET2FF::getscalarff( EvtId parent, EvtId, double t, double mass,
                              double* f0p, double* f0m )
{
    double mb = EvtPDL::getMeanMass( parent );
    double w = ( ( mb * mb ) + ( mass * mass ) - t ) / ( 2.0 * mb * mass );

    // Form factors have a general form, with parameters passed in
    // from the arguements.

    // Use dispersion relation parametrization from
    // I.Caprini, L.Lelluch, M.Neubert, Nucl. Phys. B 530,153(1998)
    const double z = ( sqrt( w + 1 ) - sqrt( 2. ) ) /
                     ( sqrt( w + 1 ) + sqrt( 2. ) );
    double v1 = v1_1 * ( 1. - 8. * rho2 * z + ( 51. * rho2 - 10. ) * z * z -
                         ( 252. * rho2 - 84. ) * z * z * z );

    *f0p = v1;
    *f0m = 0.0;

    if ( extended ) {
        //if in extended mode, use helicity-suppressed FF using the result from
        //Tanaka and Watanabe, Phys. Rev. D 82 034027 (2010)

        double r = mass / mb;
        double r1 = r + 1.0;
        double s1 = v1 * ( 1 + delta * ( -0.019 + 0.041 * ( w - 1 ) -
                                         0.015 * ( w - 1 ) *
                                             ( w - 1 ) ) );    //as in ref
        *f0m = s1 * ( w + 1.0 ) * 2.0 * r /
               ( r1 * r1 );    //convert to convention used by EvtGen with proper relative normalization
    }
}

void EvtHQET2FF::getvectorff( EvtId parent, EvtId, double t, double mass,
                              double* a1f, double* a2f, double* vf, double* a0f )
{
    double mb = EvtPDL::getMeanMass( parent );
    double w = ( ( mb * mb ) + ( mass * mass ) - t ) / ( 2.0 * mb * mass );

    // Form factors have a general form, with parameters passed in
    // from the arguements.

    double rstar = ( 2.0 * sqrt( mb * mass ) ) / ( mb + mass );

    // Use dispersion relation parametrization from
    // I.Caprini, L.Lelluch, M.Neubert, Nucl. Phys. B 530,153(1998)
    const double z = ( sqrt( w + 1 ) - sqrt( 2. ) ) /
                     ( sqrt( w + 1 ) + sqrt( 2. ) );
    double ha1 = ha1_1 * ( 1. - 8. * rho2 * z + ( 53. * rho2 - 15. ) * z * z -
                           ( 231. * rho2 - 91. ) * z * z * z );
    double r1 = r1_1 - 0.12 * ( w - 1 ) + 0.05 * ( w - 1 ) * ( w - 1 );
    double r2 = r2_1 + 0.11 * ( w - 1 ) - 0.06 * ( w - 1 ) * ( w - 1 );

    *a1f = ( 1.0 - ( t / ( ( mb + mass ) * ( mb + mass ) ) ) ) * ha1;
    *a1f = ( *a1f ) / rstar;
    *a2f = ( r2 / rstar ) * ha1;
    *vf = ( r1 / rstar ) * ha1;
    *a0f = 0.0;

    if ( extended ) {
        // Here we use the expectation for R_0(w) from
        // Fajfer et al, Phys. Rev. D 85 094025 (2012)

        double r0 = r0_1 - 0.11 * ( w - 1 ) + 0.01 * ( w - 1 ) * ( w - 1 );
        *a0f = ( r0 / rstar ) * ha1;
    }
}

void EvtHQET2FF::gettensorff( EvtId, EvtId, double, double, double*, double*,
                              double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :gettensorff in EvtHQET2FF.\n";
    ::abort();
}

void EvtHQET2FF::getbaryonff( EvtId, EvtId, double, double, double*, double*,
                              double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtHQET2FF.\n";
    ::abort();
}

void EvtHQET2FF::getdiracff( EvtId, EvtId, double, double, double*, double*,
                             double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtHQET2FF.\n";
    ::abort();
}

void EvtHQET2FF::getraritaff( EvtId, EvtId, double, double, double*, double*,
                              double*, double*, double*, double*, double*,
                              double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtHQET2FF.\n";
    ::abort();
}
