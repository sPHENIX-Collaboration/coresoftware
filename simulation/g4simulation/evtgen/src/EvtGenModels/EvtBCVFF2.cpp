
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

#include "EvtGenModels/EvtBCVFF2.hh"

using namespace std;

EvtBCVFF2::EvtBCVFF2( int idV, int fit )
{
    idVector = idV;
    whichfit = fit;
    //cout<<"==== EvtBCVFF2:: idVector="<<idVector<<" whichfit="<<whichfit<<endl;
    return;
}

void EvtBCVFF2::getvectorff( EvtId, EvtId, double t, double, double* a1f,
                             double* a2f, double* vf, double* a0f )
{
    double q2 = t;

    if ( whichfit == 0 ) {
        *vf = 0;
        *a0f = 0;
        *a1f = 1;
        *a2f = 0;

        return;
    };

    if ( idVector == EvtPDL::getId( "J/psi" ).getId() ) {    // Bc -> J/psi
        if ( whichfit == 1 ) {    // SR form factor set from [Kiselev, hep-ph/0211021]
            double Mbc = 6.277, Mpsi = 3.0967;    // Experimental values
            double Mpole2 = 4.5 * 4.5, den = 1. / ( 1. - q2 / Mpole2 );
            double FV = 0.11 * den, FAp = -0.071 * den, FA0 = 5.9 * den,
                   FAm = 0.12 * den;
            *vf = ( Mbc + Mpsi ) * FV;
            *a2f = -( Mbc + Mpsi ) * FAp;
            *a1f = FA0 / ( Mbc + Mpsi );
            *a0f = ( q2 * FAm + ( Mbc + Mpsi ) * ( *a1f ) -
                     ( Mbc - Mpsi ) * ( *a2f ) ) /
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
        };
    } else if ( idVector ==
                EvtPDL::getId( "psi(2S)" ).getId() ) {    // Bc -> psi((2S)
        if ( whichfit == 1 ) {
            ////cout<<"BC2:: psi2S, Kiselev, q2="<<q2<<endl;
            double Mbc = 6.277, Mpsi = 3.0967, Mpsi2S = 3.686,
                   kappa = Mpsi / Mpsi2S;    // Experimental values
            double Mpole2 = 4.5 * 4.5, den = 1. / ( 1. - q2 / Mpole2 );
            double FV = 0.11 * den * kappa / 3.1,
                   FAp = -0.071 * den * kappa / 4.9,
                   FA0 = 5.9 * den * kappa / 3.5, FAm = 0.12 * den * kappa / 2.3;
            *vf = ( Mbc + Mpsi2S ) * FV;
            *a2f = -( Mbc + Mpsi2S ) * FAp;
            *a1f = FA0 / ( Mbc + Mpsi2S );
            *a0f = ( q2 * FAm + ( Mbc + Mpsi2S ) * ( *a1f ) -
                     ( Mbc - Mpsi2S ) * ( *a2f ) ) /
                   ( 2 * Mpsi2S );
            return;
        } else if ( whichfit == 2 ) {
            ////cout<<"BC2:: psi2S, Ebert, q2="<<q2<<endl;
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
        };
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Not implemented :getbaryonff in EvtBCVFF2.\n";
        ::abort();
    };
}

void EvtBCVFF2::getscalarff( EvtId, EvtId, double, double, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtBCVFF2.\n";
    ::abort();
}

void EvtBCVFF2::gettensorff( EvtId, EvtId, double, double, double*, double*,
                             double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtBCVFF2.\n";
    ::abort();
}

void EvtBCVFF2::getbaryonff( EvtId, EvtId, double, double, double*, double*,
                             double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtBCVFF2.\n";
    ::abort();
}

void EvtBCVFF2::getdiracff( EvtId, EvtId, double, double, double*, double*,
                            double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtBCVFF2.\n";
    ::abort();
}

void EvtBCVFF2::getraritaff( EvtId, EvtId, double, double, double*, double*,
                             double*, double*, double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtBCVFF2.\n";
    ::abort();
}
