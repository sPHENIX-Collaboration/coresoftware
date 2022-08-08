
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

#include "EvtGenModels/EvtHQETFF.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <math.h>
#include <stdlib.h>
#include <string>

EvtHQETFF::EvtHQETFF( double hqetrho2, double hqetr1, double hqetr2,
                      double quadTerm )
{
    rho2 = hqetrho2;
    r1 = hqetr1;
    r2 = hqetr2;
    c = quadTerm;

    return;
}

EvtHQETFF::EvtHQETFF( double hqetrho2, double quadTerm )
{
    rho2 = hqetrho2;
    c = quadTerm;

    return;
}

void EvtHQETFF::getscalarff( EvtId parent, EvtId, double t, double mass,
                             double* f0p, double* f0m )
{
    double mb = EvtPDL::getMeanMass( parent );
    double w = ( ( mb * mb ) + ( mass * mass ) - t ) / ( 2.0 * mb * mass );

    // Form factors have a general form, with parameters passed in
    // from the arguements.

    double ha1 = 1 - rho2 * ( w - 1 ) + c * ( w - 1 ) * ( w - 1 );

    *f0p = ha1;
    *f0m = 0.0;

    return;
}

void EvtHQETFF::getvectorff( EvtId parent, EvtId, double t, double mass,
                             double* a1f, double* a2f, double* vf, double* a0f )
{
    double mb = EvtPDL::getMeanMass( parent );
    double w = ( ( mb * mb ) + ( mass * mass ) - t ) / ( 2.0 * mb * mass );

    // Form factors have a general form, with parameters passed in
    // from the arguements.

    double rstar = ( 2.0 * sqrt( mb * mass ) ) / ( mb + mass );
    double ha1 = 1 - rho2 * ( w - 1 );

    *a1f = ( 1.0 - ( t / ( ( mb + mass ) * ( mb + mass ) ) ) ) * ha1;
    *a1f = ( *a1f ) / rstar;
    *a2f = ( r2 / rstar ) * ha1;
    *vf = ( r1 / rstar ) * ha1;
    *a0f = 0.0;

    return;
}

void EvtHQETFF::gettensorff( EvtId, EvtId, double, double, double*, double*,
                             double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :gettensorff in EvtHQETFF.\n";
    ::abort();
}

void EvtHQETFF::getbaryonff( EvtId, EvtId, double, double, double*, double*,
                             double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtHQETFF.\n";
    ::abort();
}

void EvtHQETFF::getdiracff( EvtId, EvtId, double, double, double*, double*,
                            double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtHQETFF.\n";
    ::abort();
}

void EvtHQETFF::getraritaff( EvtId, EvtId, double, double, double*, double*,
                             double*, double*, double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtHQETFF.\n";
    ::abort();
}
