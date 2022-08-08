
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

#include "EvtGenModels/EvtBToVlnuBallFF.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <math.h>
#include <stdlib.h>
#include <string>

EvtBToVlnuBallFF::EvtBToVlnuBallFF( double r2_A1, double mfit2_A1, double r1_A2,
                                    double r2_A2, double mfit2_A2, double r1_V,
                                    double r2_V, double mfit2_V )
{
    _r2_A1 = r2_A1;
    _mfit2_A1 = mfit2_A1;
    _r1_A2 = r1_A2;
    _r2_A2 = r2_A2;
    _mfit2_A2 = mfit2_A2;
    _r1_V = r1_V;
    _r2_V = r2_V;
    _mfit2_V = mfit2_V;

    return;
}

void EvtBToVlnuBallFF::getvectorff( EvtId parent, EvtId /*daught*/, double t,
                                    double /*mass*/, double* a1f, double* a2f,
                                    double* vf, double* a0f )
{
    // FF calculations taken from the LCSR calculation of
    // P. Ball, R. Zwicky, Phys.~Rev.~{\bf D71} 014029 (2005), hep-ph/0412079.

    //Define mBstar
    EvtId Bplus = EvtPDL::getId( "B+" );
    EvtId Bminus = EvtPDL::getId( "B-" );
    double mBstar = EvtPDL::getMeanMass( EvtPDL::getId( "B*0" ) );
    if ( parent == Bplus || parent == Bminus )
        mBstar = EvtPDL::getMeanMass( EvtPDL::getId( "B*+" ) );

    double q2 = t;
    *a1f = _r2_A1 / ( 1. - q2 / _mfit2_A1 );
    *a2f = _r1_A2 / ( 1. - q2 / _mfit2_A2 ) +
           _r2_A2 / pow( 1. - q2 / _mfit2_A2, 2. );
    *vf = _r1_V / ( 1. - q2 / mBstar / mBstar ) + _r2_V / ( 1. - q2 / _mfit2_V );
    *a0f = 0.0;

    return;

    // OLD STUFF from HQETFF

    //   double mb=EvtPDL::getMeanMass(parent);
    //   double w = ((mb*mb)+(mass*mass)-t)/(2.0*mb*mass);

    // Form factors have a general form, with parameters passed in
    // from the arguements.

    //   double rstar = ( 2.0*sqrt(mb*mass))/(mb+mass);
    //   double ha1 = 1-rho2*(w-1);

    //   *a1f = (1.0 - (t/((mb+mass)*(mb+mass))))*ha1;
    //   *a1f = (*a1f)/rstar;
    //   *a2f = (r2/rstar)*ha1;
    //   *vf = (r1/rstar)*ha1;
}

void EvtBToVlnuBallFF::getscalarff( EvtId, EvtId, double, double, double*,
                                    double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getvectorff in EvtBToVlnuBallFF.\n";
    ::abort();
}

void EvtBToVlnuBallFF::gettensorff( EvtId, EvtId, double, double, double*,
                                    double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :gettensorff in EvtBToVlnuBallFF.\n";
    ::abort();
}

void EvtBToVlnuBallFF::getbaryonff( EvtId, EvtId, double, double, double*,
                                    double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtBToVlnuBallFF.\n";
    ::abort();
}

void EvtBToVlnuBallFF::getdiracff( EvtId, EvtId, double, double, double*,
                                   double*, double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtBToVlnuBallFF.\n";
    ::abort();
}

void EvtBToVlnuBallFF::getraritaff( EvtId, EvtId, double, double, double*,
                                    double*, double*, double*, double*, double*,
                                    double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtBToVlnuBallFF.\n";
    ::abort();
}
