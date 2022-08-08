
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

#include "EvtGenModels/EvtBCVFF.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

using namespace std;

EvtBCVFF::EvtBCVFF( int idV, int fit )
{
    idVector = idV;
    whichfit = fit;
    MBc = EvtPDL::getMeanMass( EvtPDL::getId( "B_c+" ) );
    MD0 = EvtPDL::getMeanMass( EvtPDL::getId( "D*0" ) );
    Mpsi = EvtPDL::getMeanMass( EvtPDL::getId( "J/psi" ) );
    Mpsi2S = EvtPDL::getMeanMass( EvtPDL::getId( "psi(2S)" ) );
    kappa = Mpsi / Mpsi2S;
    Mchi = EvtPDL::getMeanMass( EvtPDL::getId( "chi_c1" ) );
    return;
}

void EvtBCVFF::getvectorff( EvtId, EvtId, double t, double, double* a1f,
                            double* a2f, double* vf, double* a0f )
{
    double q2 = t;

    if ( whichfit == 0 ) {
        *vf = 0;
        *a0f = 0;
        *a1f = 1;
        *a2f = 0;
        return;
    }

    if ( idVector == EvtPDL::getId( "J/psi" ).getId() ) {    // Bc -> J/psi
        if ( whichfit == 1 ) {    // SR form factor set from [Kiselev, hep-ph/0211021]
            double Mpole2 = 4.5 * 4.5, den = 1. / ( 1. - q2 / Mpole2 );
            double FV = 0.11 * den, FAp = -0.074 * den, FA0 = 5.9 * den,
                   FAm = 0.12 * den;
            *vf = ( MBc + Mpsi ) * FV;
            *a2f = -( MBc + Mpsi ) * FAp;
            *a1f = FA0 / ( MBc + Mpsi );
            *a0f = ( q2 * FAm + ( MBc + Mpsi ) * ( *a1f ) -
                     ( MBc - Mpsi ) * ( *a2f ) ) /
                   ( 2 * Mpsi );
            return;
        } else if ( whichfit ==
                    2 ) {    // form factor set from  [Ebert, hep-ph/0306306]
            *vf = ( 0.49077824756158533 - 0.0012925655191347828 * q2 ) /
                  ( 1 - 0.06292520325875656 * q2 );
            *a0f = ( 0.4160345034630221 - 0.0024720095310225023 * q2 ) /
                   ( 1 - 0.061603451915567785 * q2 );
            *a1f = ( 0.4970212860605933 - 0.0067519730024654745 * q2 ) /
                   ( 1 - 0.050487026667172176 * q2 );
            *a2f = ( 0.7315284919705497 + 0.0014263826220727142 * q2 -
                     0.0006946090066269195 * q2 * q2 ) /
                   ( 1 - 0.04885587273651653 * q2 );
            return;
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtBCVFF" )
                << "Must choose 0 (a1f = 1), 1 (Kiselev), or 2 (Ebert).\n";
            ::abort();
        }
    } else if ( idVector ==
                EvtPDL::getId( "psi(2S)" ).getId() ) {    // Bc -> psi((2S)
        if ( whichfit == 1 ) {
            double Mpole2 = 4.5 * 4.5, den = 1. / ( 1. - q2 / Mpole2 );
            double FV = 0.11 * den * kappa / 3.1,
                   FAp = -0.074 * den * kappa / 4.9,
                   FA0 = 5.9 * den * kappa / 3.5, FAm = 0.12 * den * kappa / 2.3;
            *vf = ( MBc + Mpsi2S ) * FV;
            *a2f = -( MBc + Mpsi2S ) * FAp;
            *a1f = FA0 / ( MBc + Mpsi2S );
            *a0f = ( q2 * FAm + ( MBc + Mpsi2S ) * ( *a1f ) -
                     ( MBc - Mpsi2S ) * ( *a2f ) ) /
                   ( 2 * Mpsi2S );
            return;
        } else if ( whichfit == 2 ) {
            *vf = ( 0.24177223968739653 - 0.053589051007278135 * q2 ) /
                  ( 1 - 0.0977848994260899 * q2 );
            *a0f = ( 0.23996026570086615 - 0.03530198514007337 * q2 ) /
                   ( 1 - 0.09371162519983989 * q2 );
            *a1f = ( 0.17418379258849329 - 0.004129699022085851 * q2 * q2 ) /
                   ( 1 + 0.06607665248402918 * q2 );
            *a2f = ( 0.1352376939112041 - 0.040361722565209444 * q2 +
                     0.003343515369431853 * q2 * q2 ) /
                   ( 1 - 0.1463698128333418 * q2 );
            return;
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtBCVFF" )
                << "Must choose 0 (a1f = 1), 1 (Kiselev), or 2 (Ebert).\n";
            ::abort();
        }
    } else if ( idVector == EvtPDL::getId( "chi_c1" ).getId() ) {    // Bc -> chi_c1
        if ( whichfit == 3 ) {    // FF from Wang et al 10.1103/PhysRevD.79.114018
            double SoverD = ( MBc + Mchi ) / ( MBc - Mchi );
            double DoverS = ( MBc - Mchi ) / ( MBc + Mchi );
            double ratio = q2 / ( MBc * MBc );

            double vf_0 = SoverD * 0.36;
            double vf_c1 = 1.98;
            double vf_c2 = 0.43;

            double a2f_0 = SoverD * 0.15;
            double a2f_c1 = 1.22;
            double a2f_c2 = -0.08;

            double a1f_0 = DoverS * 0.85;
            double a1f_c1 = -0.51;
            double a1f_c2 = -1.38;

            double a0f_0 = 0.13;
            double a0f_c1 = 2.99;
            double a0f_c2 = 0.023;

            *vf = vf_0 * exp( vf_c1 * ratio + vf_c2 * ratio * ratio );
            *a2f = a2f_0 * exp( a2f_c1 * ratio + a2f_c2 * ratio * ratio );
            *a1f = a1f_0 * exp( a1f_c1 * ratio + a1f_c2 * ratio * ratio );
            *a0f = a0f_0 * exp( a0f_c1 * ratio + a0f_c2 * ratio * ratio );
            return;
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtBCVFF" )
                << "Must choose 0 (a1f = 1) or 3 (Wang).\n";
            ::abort();
        }
    } else if ( idVector == EvtPDL::getId( "D*0" ).getId() ||
                idVector == EvtPDL::getId( "anti-D*0" ).getId() ) {
        if ( whichfit == 1 ) {
            // SR form factor set from Kiselev, hep-ph/0211021
            double Mpole2 = 6.2 * 6.2, den = ( 1. - q2 / Mpole2 );
            if ( fabs( den ) < 1e-10 ) {
                *vf = 0;
                *a2f = 0;
                *a1f = 0;
                *a0f = 0;
            } else {
                double FV = 0.20 / den, FAp = -0.062 / den, FA0 = 3.6,
                       FAm = 0.11 / den;
                *vf = ( MBc + MD0 ) * FV;
                *a2f = -( MBc + MD0 ) * FAp;
                *a1f = FA0 / ( MBc + MD0 );
                *a0f = ( q2 * FAm + ( MBc + MD0 ) * ( *a1f ) -
                         ( MBc - MD0 ) * ( *a2f ) ) /
                       ( 2 * MD0 );
            }
            return;
        } else if ( whichfit == 2 ) {
            // form factors from Ebert, hep-ph/0306306
            double ratio = q2 / MBc / MBc;
            double const fV_0 = 0.202, fV_a = 1.38, fV_b = 1.31;
            double const fA2_0 = 0.22, fA2_a = 2.44, fA2_b = -1.21;
            double const fA0_0 = 0.144, fA0_a = 1.18, fA0_b = 1.39;
            double const fA1_0 = 0.174, fA1_a = 1.69, fA1_b = -0.219;

            *vf = fV_0 / ( 1 - fV_a * ratio - fV_b * ratio * ratio );
            *a2f = fA2_0 / ( 1 - fA2_a * ratio - fA2_b * ratio * ratio );
            *a0f = fA0_0 / ( 1 - fA0_a * ratio - fA0_b * ratio * ratio );
            *a1f = fA1_0 / ( 1 - fA1_a * ratio - fA1_b * ratio * ratio );
            return;
        } else {
            EvtGenReport( EVTGEN_ERROR, "BCVFF" )
                << "FF should be 0 (a1f = 1), 1 (Kiselev 2002) or 2 (Ebert).\n";
            ::abort();
        }
    } else {
        EvtGenReport( EVTGEN_ERROR, "ECVFF" )
            << "Only J/psi, psi(2S), chi_c1 and D*0/anti-D*0 implemented.\n";
        ::abort();
    }
}

void EvtBCVFF::getscalarff( EvtId, EvtId, double, double, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getscalarff in EvtBCVFF.\n";
    ::abort();
}

void EvtBCVFF::gettensorff( EvtId, EvtId, double, double, double*, double*,
                            double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :gettensorff in EvtBCVFF.\n";
    ::abort();
}

void EvtBCVFF::getbaryonff( EvtId, EvtId, double, double, double*, double*,
                            double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtBCVFF.\n";
    ::abort();
}

void EvtBCVFF::getdiracff( EvtId, EvtId, double, double, double*, double*,
                           double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtBCVFF.\n";
    ::abort();
}

void EvtBCVFF::getraritaff( EvtId, EvtId, double, double, double*, double*,
                            double*, double*, double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtBCVFF.\n";
    ::abort();
}
