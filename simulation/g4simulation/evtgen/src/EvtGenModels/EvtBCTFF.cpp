
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

#include "EvtGenModels/EvtBCTFF.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

using namespace std;

EvtBCTFF::EvtBCTFF( int idT, int fit )
{
    idTensor = idT;
    whichfit = fit;
    MBc = EvtPDL::getMeanMass( EvtPDL::getId( "B_c+" ) );
    return;
}

void EvtBCTFF::getscalarff( EvtId, EvtId, double, double, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getscalarff in EvtBCTFF.\n";
    ::abort();
}

void EvtBCTFF::getvectorff( EvtId, EvtId, double, double, double*, double*,
                            double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getvectorff in EvtBCTFF.\n";
    ::abort();
}

void EvtBCTFF::gettensorff( EvtId /*p*/, EvtId /*d*/, double t, double /*mass*/,
                            double* hf, double* kf, double* bpf, double* bmf )
{
    double q2 = t;

    if ( whichfit == 0 ) {
        *hf = 0;
        *kf = 0;
        *bpf = 0;
        *bmf = 0;
        return;
    }

    if ( idTensor == EvtPDL::getId( "chi_c2" ).getId() ) {    // Bc -> chi_c1
        if ( whichfit == 3 ) {    // FF from Wang et al 10.1103/PhysRevD.79.114018
            double ratio = q2 / ( MBc * MBc );

            double hf_0 = 0.022;
            double hf_c1 = 2.58;
            double hf_c2 = 0.61;

            double kf_0 = 1.27;
            double kf_c1 = 1.61;
            double kf_c2 = 0.24;

            double bpf_0 = -0.011;
            double bpf_c1 = 2.27;
            double bpf_c2 = 0.46;

            double bmf_0 = 0.020;
            double bmf_c1 = 2.48;
            double bmf_c2 = 0.56;

            *hf = hf_0 * exp( hf_c1 * ratio + hf_c2 * ratio * ratio );
            *kf = kf_0 * exp( kf_c1 * ratio + kf_c2 * ratio * ratio );
            *bpf = bpf_0 * exp( bpf_c1 * ratio + bpf_c2 * ratio * ratio );
            *bmf = bmf_0 * exp( bmf_c1 * ratio + bmf_c2 * ratio * ratio );
            return;
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Must choose 0 (a1f = 1) or 3 (Wang).\n";
            ::abort();
        }
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "chi_c2 is the only (pseudo)vector decay implemented in EvtBCTFF.\n";
        ::abort();
    }
}

void EvtBCTFF::getbaryonff( EvtId, EvtId, double, double, double*, double*,
                            double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtBCTFF.\n";
    ::abort();
}

void EvtBCTFF::getdiracff( EvtId, EvtId, double, double, double*, double*,
                           double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtBCTFF.\n";
    ::abort();
}

void EvtBCTFF::getraritaff( EvtId, EvtId, double, double, double*, double*,
                            double*, double*, double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtBCTFF.\n";
    ::abort();
}
