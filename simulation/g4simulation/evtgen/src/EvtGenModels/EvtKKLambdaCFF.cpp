
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

#include "EvtGenModels/EvtKKLambdaCFF.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <math.h>
#include <stdlib.h>
#include <string>

EvtKKLambdaCFF::EvtKKLambdaCFF( int numarg, double* arglist )
{
    _nargs = numarg;
    for ( int i = 0; i < numarg; i++ ) {
        _args[i] = arglist[i];
    }

    return;
}

void EvtKKLambdaCFF::getbaryonff( EvtId /*parent*/, EvtId /*daught*/, double t,
                                  double /*mass*/, double* f1v, double* f1a,
                                  double* f2v, double* f2a )
{
    *f1v = ( _args[0] ) / ( 1.0 - ( t / ( _args[1] * _args[1] ) ) );

    *f2v = 0.;
    *f2a = 0.;
    *f1a = -1.0 * ( *f1v );
}

void EvtKKLambdaCFF::getscalarff( EvtId, EvtId, double, double, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getscalarff in EvtKKLambdaCFF.\n";
    ::abort();
}

void EvtKKLambdaCFF::getvectorff( EvtId, EvtId, double, double, double*,
                                  double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getvectorff in EvtKKLambdaCFF.\n";
    ::abort();
}

void EvtKKLambdaCFF::gettensorff( EvtId, EvtId, double, double, double*,
                                  double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :gettensorff in EvtKKLambdaCFF.\n";
    ::abort();
}

void EvtKKLambdaCFF::getdiracff( EvtId, EvtId, double, double, double*, double*,
                                 double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtKKLambdaCFF.\n";
    ::abort();
}

void EvtKKLambdaCFF::getraritaff( EvtId, EvtId, double, double, double*,
                                  double*, double*, double*, double*, double*,
                                  double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtKKLambdaCFF.\n";
    ::abort();
}
