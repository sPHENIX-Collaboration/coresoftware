
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

#include "EvtGenModels/EvtLb2BaryonlnuFF.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <math.h>
#include <stdlib.h>
#include <string>
using std::endl;

void EvtLb2BaryonlnuFF::getdiracff( EvtId parent, EvtId daught, double q2,
                                    double /* mass */, double* f1, double* f2,
                                    double* f3, double* g1, double* g2,
                                    double* g3 )
{
    // Define Event IDs for Lb and p, N+ and Lc+ states
    static EvtId LAMB = EvtPDL::getId( "Lambda_b0" );
    static EvtId LAMBB = EvtPDL::getId( "anti-Lambda_b0" );
    static EvtId PRO = EvtPDL::getId( "p+" );
    static EvtId PROB = EvtPDL::getId( "anti-p-" );
    static EvtId N1440 = EvtPDL::getId( "N(1440)+" );
    static EvtId N1440B = EvtPDL::getId( "anti-N(1440)-" );
    static EvtId N1535 = EvtPDL::getId( "N(1535)+" );
    static EvtId N1535B = EvtPDL::getId( "anti-N(1535)-" );
    static EvtId N1650 = EvtPDL::getId( "N(1650)+" );
    static EvtId N1650B = EvtPDL::getId( "anti-N(1650)-" );
    static EvtId N1710 = EvtPDL::getId( "N(1710)+" );
    static EvtId N1710B = EvtPDL::getId( "anti-N(1710)-" );
    static EvtId LAMCP = EvtPDL::getId( "Lambda_c+" );
    static EvtId LAMCM = EvtPDL::getId( "anti-Lambda_c-" );
    static EvtId LAMC1P = EvtPDL::getId( "Lambda_c(2593)+" );
    static EvtId LAMC1M = EvtPDL::getId( "anti-Lambda_c(2593)-" );

    double F1, F2, F3, G1, G2, G3;

    if ( ( parent == LAMB && daught == PRO ) ||
         ( parent == LAMBB && daught == PROB ) ||
         ( parent == LAMB && daught == LAMCP ) ||
         ( parent == LAMBB && daught == LAMCM ) ) {
        //  Parameters needed in the calculation;
        double mQ = 5.28;
        double aL = 0.59;
        double md = 0.40;
        double MLamB = EvtPDL::getMass( parent );
        double MLamq = EvtPDL::getMass( daught );
        double mq = md;
        double aLp = 0.48;

        //set mq and aLp based on whether Lb->Lc* or Lb->N*
        if ( ( parent == LAMB && daught == LAMCP ) ||
             ( parent == LAMBB && daught == LAMCM ) ) {
            mq = 1.89;
            aLp = 0.55;
        }

        double aL2 = aL * aL;
        double aLp2 = aLp * aLp;
        double aLLp2 = 0.5 * ( aL2 + aLp2 );

        // relativistic correction factor
        double k2 = 1.0;
        double rho2 = 3. * md * md / ( 2. * k2 * aLLp2 );

        // w = scalar product of the 4 velocities of the Lb and Lc.
        double w = 0.5 * ( MLamB * MLamB + MLamq * MLamq - q2 ) / MLamB / MLamq;

        double I = pow( aL * aLp / aLLp2, 1.5 ) * exp( -rho2 * ( w * w - 1. ) );

        // Calculate the form factors
        F1 = I * ( 1.0 + ( md / aLLp2 ) * ( ( aLp2 / mq ) + ( aL2 / mQ ) ) );
        F2 = -I * ( ( md / mq ) * ( aLp2 / aLLp2 ) -
                    aL2 * aLp2 / ( 4. * aLLp2 * mq * mQ ) );
        F3 = -I * md * aL2 / ( mQ * aLLp2 );

        G1 = I * ( 1.0 - ( aL2 * aLp2 ) / ( 12. * aLLp2 * mq * mQ ) );
        G2 = -I * ( md * aLp2 / ( mq * aLLp2 ) +
                    ( aL2 * aLp2 ) / ( 12. * aLLp2 * mq * mQ ) *
                        ( 1. + 12. * md * md / aLLp2 ) );
        G3 = I * ( md * aL2 / ( mQ * aLLp2 ) +
                   md * md * aL2 * aLp2 / ( mq * mQ * aLLp2 * aLLp2 ) );

        // Set form factors to be passed to the amplitude calc.
        *f1 = F1;
        *f2 = F2;
        *f3 = F3;
        *g1 = G1;
        *g2 = G2;
        *g3 = G3;
    } else if ( ( parent == LAMB && daught == N1440 ) ||
                ( parent == LAMBB && daught == N1440B ) ||
                ( parent == LAMB && daught == N1710 ) ||
                ( parent == LAMBB && daught == N1710B ) ) {
        //  Parameters needed in the calculation;
        double mQ = 5.28;
        double md = 0.40;
        double mq = md;
        double MLamB = EvtPDL::getMass( parent );
        double MLamq = EvtPDL::getMass( daught );

        double aL = 0.59;
        double aLp = 0.48;

        double aL2 = aL * aL;
        double aLp2 = aLp * aLp;
        double aLLp2 = 0.5 * ( aL2 + aLp2 );

        // relativistic correction factor
        double k2 = 1.0;
        double rho2 = 3. * md * md / ( 2. * k2 * aLLp2 );

        // w = scalar product of the 4 velocities of the Lb and Lc.
        double w = 0.5 * ( MLamB * MLamB + MLamq * MLamq - q2 ) / MLamB / MLamq;

        double I = sqrt( 1.5 ) * pow( aL * aLp / aLLp2, 1.5 ) *
                   exp( -rho2 * ( w * w - 1. ) );

        // Calculate the form factors
        F1 = ( I / ( 2. * aLLp2 ) ) *
             ( ( aL2 - aLp2 ) - ( md / ( 3. * aLLp2 ) ) *
                                    ( ( aLp2 / mq ) * ( 7. * aL2 - 3. * aLp2 ) +
                                      ( aL2 / mQ ) * ( 7. * aLp2 - 3. * aL2 ) ) );

        F2 = -I * ( aLp2 / ( 6. * mq * aLLp2 * aLLp2 ) ) *
             ( 7. * aL2 - 3. * aLp2 ) * ( md - ( aL2 / ( 4. * mQ ) ) );

        F3 = I * ( aL2 * md / ( 6. * mQ * aLLp2 * aLLp2 ) ) *
             ( 7. * aLp2 - 3. * aL2 );

        G1 = I * ( ( aL2 - aLp2 ) / ( 2 * aLLp2 ) -
                   ( aL2 * aLp2 * ( 7. * aL2 - 3. * aLp2 ) ) /
                       ( 72. * aLLp2 * aLLp2 * mq * mQ ) );

        G2 = -I * ( aLp2 / ( 6. * mq * aLLp2 * aLLp2 ) ) *
             ( ( 7. * aL2 - 3. * aLp2 ) * ( md + ( aL2 / ( 6. * mQ ) ) ) +
               ( 7. * md * md * aL2 * ( aL2 - aLp2 ) / ( mQ * aLLp2 ) ) );

        G3 = -I * ( aL2 * md / ( 6. * mQ * aLLp2 * aLLp2 ) ) *
             ( ( 7. * aLp2 - 3. * aL2 ) -
               ( 7 * md * aLp2 * ( aL2 - aLp2 ) / ( mq * aLLp2 ) ) );

        // Set form factors to be passed to the amplitude calc.
        *f1 = F1;
        *f2 = F2;
        *f3 = F3;
        *g1 = G1;
        *g2 = G2;
        *g3 = G3;

    } else if ( ( parent == LAMB && daught == N1535 ) ||
                ( parent == LAMBB && daught == N1535B ) ||
                ( parent == LAMB && daught == N1650 ) ||
                ( parent == LAMBB && daught == N1650B ) ||
                ( parent == LAMB && daught == LAMC1P ) ||
                ( parent == LAMBB && daught == LAMC1M ) ) {
        double mQ = 5.28;
        double md = 0.40;
        double aL = 0.59;
        double MLamB = EvtPDL::getMass( parent );
        double MLamq = EvtPDL::getMass( daught );
        double mq = md;
        double aLp = 0.37;

        //set mq and aLp based on whether Lb->Lc* or Lb->N*
        if ( ( parent == LAMB && daught == LAMC1P ) ||
             ( parent == LAMBB && daught == LAMC1M ) ) {
            mq = 1.89;
            aLp = 0.47;
        }

        double aL2 = aL * aL;
        double aLp2 = aLp * aLp;
        double aLLp2 = 0.5 * ( aL2 + aLp2 );

        // relativistic correction factor
        double k2 = 1.0;
        double rho2 = 3. * md * md / ( 2. * k2 * aLLp2 );

        // w = scalar product of the 4 velocities of the Lb and Lc.
        double w = 0.5 * ( MLamB * MLamB + MLamq * MLamq - q2 ) / MLamB / MLamq;

        double I = pow( aL * aLp / aLLp2, 2.5 ) * exp( -rho2 * ( w * w - 1. ) );

        // Calculate the form factors
        F1 = I * aL / 6.0 * ( 3.0 / mq - 1.0 / mQ );
        F2 = -I * ( 2.0 * md / aL - aL / ( 2.0 * mq ) +
                    2. * md * md * aL / ( mQ * aLLp2 ) -
                    ( md * aL / ( 6. * mq * mQ * aLLp2 ) ) *
                        ( 3. * aL2 - 2. * aLp2 ) );
        F3 = I * 2. * md * md * aL / ( mQ * aLLp2 );

        G1 = I * ( 2.0 * md / aL - aL / ( 6. * mQ ) +
                   ( md * aL / ( 6. * mq * mQ * aLLp2 ) ) *
                       ( 3. * aL2 - 2. * aLp2 ) );
        G2 = I * ( -2. * md / aL + aL / ( 2. * mq ) + aL / ( 3. * mQ ) );
        G3 = I * aL / ( 3. * mQ ) *
             ( 1.0 - ( md / ( 2. * mq * aLLp2 ) ) * ( 3. * aL2 - 2. * aLp2 ) );

        // Set form factors to be passed to the amplitude calc.
        *f1 = F1;
        *f2 = F2;
        *f3 = F3;
        *g1 = G1;
        *g2 = G2;
        *g3 = G3;
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Only Lb -> N*+ transitions allowed in EvtLb2BaryonlnuFF.\n";
        ::abort();
    }

    return;
}

void EvtLb2BaryonlnuFF::getraritaff( EvtId parent, EvtId daught, double q2,
                                     double, double* f1, double* f2, double* f3,
                                     double* f4, double* g1, double* g2,
                                     double* g3, double* g4 )
{
    static EvtId LAMB = EvtPDL::getId( "Lambda_b0" );
    static EvtId LAMBB = EvtPDL::getId( "anti-Lambda_b0" );
    static EvtId N1520 = EvtPDL::getId( "N(1520)+" );
    static EvtId N1520B = EvtPDL::getId( "anti-N(1520)-" );
    static EvtId N1720 = EvtPDL::getId( "N(1720)+" );
    static EvtId N1720B = EvtPDL::getId( "anti-N(1720)-" );
    static EvtId N1700 = EvtPDL::getId( "N(1700)+" );
    static EvtId N1700B = EvtPDL::getId( "anti-N(1700)-" );
    static EvtId N1900 = EvtPDL::getId( "N(1900)+" );
    static EvtId N1900B = EvtPDL::getId( "anti-N(1900)-" );
    static EvtId N1875 = EvtPDL::getId( "N(1875)+" );
    static EvtId N1875B = EvtPDL::getId( "anti-N(1875)-" );
    static EvtId LAMC2P = EvtPDL::getId( "Lambda_c(2625)+" );
    static EvtId LAMC2M = EvtPDL::getId( "anti-Lambda_c(2625)-" );

    double F1, F2, F3, F4, G1, G2, G3, G4;

    // 3/2 - case
    if ( ( parent == LAMB && daught == N1520 ) ||
         ( parent == LAMBB && daught == N1520B ) ||
         ( parent == LAMB && daught == N1700 ) ||
         ( parent == LAMBB && daught == N1700B ) ||
         ( parent == LAMB && daught == N1875 ) ||
         ( parent == LAMBB && daught == N1875B ) ||
         ( parent == LAMB && daught == LAMC2P ) ||
         ( parent == LAMBB && daught == LAMC2M ) ) {
        double mQ = 5.28;
        double md = 0.40;
        double aL = 0.59;
        double MLamB = EvtPDL::getMass( parent );
        double MLamq = EvtPDL::getMass( daught );
        double mq = md;
        double aLp = 0.37;

        //set mq and aLp based on whether Lb->Lc* or Lb->N*
        if ( ( parent == LAMB && daught == LAMC2P ) ||
             ( parent == LAMBB && daught == LAMC2M ) ) {
            mq = 1.89;
            aLp = 0.47;
        }

        double aL2 = aL * aL;
        double aLp2 = aLp * aLp;
        double aLLp2 = 0.5 * ( aL2 + aLp2 );

        // relativistic correction factor
        double k2 = 1.0;
        double rho2 = 3. * md * md / ( 2. * k2 * aLLp2 );

        // w = scalar product of the 4 velocities of the Lb and Lc.
        double w = 0.5 * ( MLamB * MLamB + MLamq * MLamq - q2 ) / MLamB / MLamq;

        double I = -( 1. / sqrt( 3. ) ) * pow( aL * aLp / aLLp2, 2.5 ) *
                   exp( -rho2 * ( w * w - 1. ) );

        // Calculate the form factors
        F1 = I * 3.0 * md / aL *
             ( 1.0 + ( md / aLLp2 ) * ( ( aLp2 / mq ) + ( aL2 / mQ ) ) );
        F2 = -I * ( ( 3. * md * md / mq ) * ( aLp2 / ( aLLp2 * aL2 ) ) -
                    5. * aL * aLp2 * md / ( 4. * aLLp2 * mq * mQ ) );
        F3 = -I * ( 3. * md * md * aL / ( mQ * aLLp2 ) + aL / ( 2. * mQ ) );
        F4 = I * aL / mQ;

        G1 = I * ( 3.0 * md / aL -
                   ( aL / ( 2. * mQ ) ) *
                       ( 1. + 3. * md * aLp2 / ( 2. * aLLp2 * mq ) ) );
        G2 = -I * ( ( 3. * md * md / mq ) * ( aLp2 / ( aLLp2 * aL ) ) +
                    aL * aLp2 * md / ( 4. * aLLp2 * aLLp2 * mq * mQ ) *
                        ( aLLp2 + 12. * md * md ) );
        G3 = I * aL / ( mQ * aLLp2 ) *
             ( aLLp2 / 2. + 3. * md * md +
               aLp2 * md / ( mq * aLLp2 ) * ( aLLp2 + 6. * md * md ) );
        G4 = -I * ( aL / mQ + md / ( mq * mQ ) * aLp2 * aL / aLLp2 );

        // Set form factors to be passed to the amplitude calc.
        *f1 = F1;
        *f2 = F2;
        *f3 = F3;
        *f4 = F4;
        *g1 = G1;
        *g2 = G2;
        *g3 = G3;
        *g4 = G4;
    }
    // 3/2 + case
    else if ( ( parent == LAMB && daught == N1720 ) ||
              ( parent == LAMBB && daught == N1720B ) ||
              ( parent == LAMB && daught == N1900 ) ||
              ( parent == LAMBB && daught == N1900B )

    ) {
        double mQ = 5.28;
        double md = 0.40;
        double mq = md;
        double MLamB = EvtPDL::getMass( parent );
        double MLamq = EvtPDL::getMass( daught );

        double aL = 0.59;
        double aLp = 0.35;

        double aL2 = aL * aL;
        double aLp2 = aLp * aLp;
        double aLLp2 = 0.5 * ( aL2 + aLp2 );

        // relativistic correction factor
        double k2 = 1.0;
        double rho2 = 3. * md * md / ( 2. * k2 * aLLp2 );

        // w = scalar product of the 4 velocities of the Lb and Lc.
        double w = 0.5 * ( MLamB * MLamB + MLamq * MLamq - q2 ) / MLamB / MLamq;

        double I = ( 1. / sqrt( 5. ) ) * pow( aL * aLp / aLLp2, 3.5 ) *
                   exp( -rho2 * ( w * w - 1. ) );

        // Calculate the form factors
        F1 = -I * ( md / 2. ) * ( ( 5. / mq ) - ( 3. / mQ ) );

        F2 = I * ( md / aL ) *
             ( ( 6. * md / aL ) - ( 5 * aL / ( 2. * mq ) ) +
               ( 6. * md * md * aL ) / ( aLLp2 * mQ ) -
               ( md * aL * ( aL2 - 2. * aLp2 ) ) / ( 2 * aLLp2 * mq * mQ ) );

        F3 = -I * ( md / mQ ) * ( 1 + ( 6. * md * md ) / aLLp2 );

        F4 = I * ( 2. * md / mQ );

        G1 = -I *
             ( ( 6. * md * md / aL2 ) - md / ( 2. * mQ ) +
               md * md * ( 11. * aL2 - 6. * aLp2 ) / ( 6. * aLLp2 * mq * mQ ) );

        G2 = I * ( 6 * md * md / aL2 - 5 * md / ( 2.0 * mq ) - ( 2 * md ) / mQ +
                   5 * aL2 / ( 12.0 * mq * mQ ) -
                   ( 2 * md * md * aL2 ) / ( 3.0 * aLLp2 * mq * mQ ) );

        G3 = -I *
             ( ( md / ( 2. * mQ ) ) - 5 * aL2 / ( 24.0 * mq * mQ ) -
               md * md * ( 5 * aL2 - 2 * aLp2 ) / ( 4.0 * mq * mQ * aLLp2 ) );

        G4 = -I * 5. * aL2 / ( 6. * mq * mQ );

        // Set form factors to be passed to the amplitude calc.
        *f1 = F1;
        *f2 = F2;
        *f3 = F3;
        *f4 = F4;
        *g1 = G1;
        *g2 = G2;
        *g3 = G3;
        *g4 = G4;
    }

    else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Only Lb -> N*+ transitions allowed in EvtLb2BaryonlnuFF.\n";
        ::abort();
    }

    return;
}

void EvtLb2BaryonlnuFF::getscalarff( EvtId, EvtId, double, double, double*,
                                     double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getscalarff in EvtLb2BaryonlnuFF.\n";
    ::abort();
}

void EvtLb2BaryonlnuFF::getvectorff( EvtId, EvtId, double, double, double*,
                                     double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getvectorff in EvtLb2BaryonlnuFF.\n";
    ::abort();
}

void EvtLb2BaryonlnuFF::gettensorff( EvtId, EvtId, double, double, double*,
                                     double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :gettensorff in EvtLb2BaryonlnuFF.\n";
    ::abort();
}

void EvtLb2BaryonlnuFF::getbaryonff( EvtId, EvtId, double, double, double*,
                                     double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtLb2BaryonlnuFF.\n";
    ::abort();
}
