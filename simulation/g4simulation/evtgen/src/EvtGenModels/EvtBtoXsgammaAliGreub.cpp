
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

#include "EvtGenModels/EvtBtoXsgammaAliGreub.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>
using std::endl;

void EvtBtoXsgammaAliGreub::init( int nArg, double* /*args*/ )
{
    if ( ( nArg - 1 ) != 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBtoXsgamma generator model "
            << "EvtBtoXsgammaAliGreub expected "
            << "zero arguments but found: " << nArg - 1 << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

double EvtBtoXsgammaAliGreub::GetMass( int Xscode )
{
    // The special lineshape for strange hadrons X_s in b -> s gamma:
    // An 18 parameter function fitted to the theoretical mass spectrum
    // of Ali & Greub for a B meson mass of 5.279 GeV; top quark mass of
    // 174.3 GeV; strange quark mass of 0.48 GeV (tuned to give minimum
    // M_Xs of 0.64 GeV) and Fermi momentum of 265 MeV for spectator quark
    // mass of 150 MeV (from CLEO fit). Truncated at max on high side
    // and min (just above K pi or KK thresold) on low side.
    double min = 0.64;
    double max = 4.5;
    double xbox, ybox, alifit;
    double mass = 0.0;

    double par[18];
    if ( ( Xscode == 30343 ) || ( Xscode == -30343 ) || ( Xscode == 30353 ) ||
         ( Xscode == -30353 ) ) {    // Xsu or Xsd
        min = 0.6373;                //  Just above K pi threshold for Xsd/u
        //min=0.6333; //  K pi threshold for neutral Xsd
        par[0] = -2057.2380371094;
        par[1] = 2502.2556152344;
        par[2] = 1151.5632324219;
        par[3] = 0.82431584596634;
        par[4] = -4110.5234375000;
        par[5] = 8445.6757812500;
        par[6] = -3034.1894531250;
        par[7] = 1.1557708978653;
        par[8] = 1765.9311523438;
        par[9] = 1.3730158805847;
        par[10] = 0.51371538639069;
        par[11] = 2.0056934356689;
        par[12] = 37144.097656250;
        par[13] = -50296.781250000;
        par[14] = 27319.095703125;
        par[15] = -7408.0678710938;
        par[16] = 1000.8093261719;
        par[17] = -53.834449768066;
    } else if ( ( Xscode == 30363 ) || ( Xscode == -30363 ) ) {
        min = 0.9964;    // Just above KK threshold for Xss
        par[0] = -32263.908203125;
        par[1] = 57186.589843750;
        par[2] = -24230.728515625;
        par[3] = 1.1155973672867;
        par[4] = -12161.131835938;
        par[5] = 20162.146484375;
        par[6] = -7198.8564453125;
        par[7] = 1.3783323764801;
        par[8] = 1995.1691894531;
        par[9] = 1.4655895233154;
        par[10] = 0.48869228363037;
        par[11] = 2.1038570404053;
        par[12] = 55100.058593750;
        par[13] = -75201.703125000;
        par[14] = 41096.066406250;
        par[15] = -11205.986328125;
        par[16] = 1522.4024658203;
        par[17] = -82.379623413086;
    } else {
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
            << "In EvtBtoXsgammaAliGreub: Particle with id " << Xscode
            << " is not a Xss particle" << endl;
        return 0;
    }

    double boxheight = par[8];
    double boxwidth = max - min;

    while ( ( mass > max ) || ( mass < min ) ) {
        xbox = EvtRandom::Flat( boxwidth ) + min;
        ybox = EvtRandom::Flat( boxheight );
        if ( xbox < par[3] ) {
            alifit = par[0] + par[1] * xbox + par[2] * pow( xbox, 2 );
        } else if ( xbox < par[7] ) {
            alifit = par[4] + par[5] * xbox + par[6] * pow( xbox, 2 );
        } else if ( xbox < par[11] ) {
            alifit = par[8] * exp( -0.5 * pow( ( xbox - par[9] ) / par[10], 2 ) );
        } else {
            alifit = par[12] + par[13] * xbox + par[14] * pow( xbox, 2 ) +
                     par[15] * pow( xbox, 3 ) + par[16] * pow( xbox, 4 ) +
                     par[17] * pow( xbox, 5 );
        }
        if ( ybox > alifit ) {
            mass = 0.0;
        } else {
            mass = xbox;
        }
    }
    return mass;
}
