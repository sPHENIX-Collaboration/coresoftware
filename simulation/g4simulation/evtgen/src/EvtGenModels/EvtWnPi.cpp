
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

#include "EvtGenModels/EvtWnPi.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParser.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenModels/EvtTauHadnu.hh"

#include <ctype.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace std;

// W+ -> pi_ current
EvtVector4C EvtWnPi::WCurrent( EvtVector4R q1 )
{
    return q1;
}

// W+ -> pi+ pi0 current
EvtVector4C EvtWnPi::WCurrent( EvtVector4R q1, EvtVector4R q2 )
{
    return BWr( q1 + q2 ) * ( q1 - q2 );
}

// W+ -> pi+ pi+ pi- current
EvtVector4C EvtWnPi::WCurrent( EvtVector4R q1, EvtVector4R q2, EvtVector4R q3 )
{
    EvtVector4R Q = q1 + q2 + q3;
    double Q2 = Q.mass2();
    return BWa( Q ) *
           ( ( q1 - q3 ) - ( Q * ( Q * ( q1 - q3 ) ) / Q2 ) * BWr( q2 + q3 ) +
             ( q2 - q3 ) - ( Q * ( Q * ( q2 - q3 ) ) / Q2 ) * BWr( q1 + q3 ) );
}

// W+ -> pi+ pi+ pi- pi- pi+ current with symmetrization
EvtVector4C EvtWnPi::WCurrent( EvtVector4R q1, EvtVector4R q2, EvtVector4R q3,
                               EvtVector4R q4, EvtVector4R q5 )
{
    //  double Q2 = Qtot*Qtot;
    //  return q1-Qtot*(q1*Qtot)/Q2;
    EvtVector4C V = JB( q1, q2, q3, q4, q5 ) + JB( q5, q2, q3, q4, q1 ) +
                    JB( q1, q5, q3, q4, q2 ) + JB( q1, q2, q4, q3, q5 ) +
                    JB( q5, q2, q4, q3, q1 ) + JB( q1, q5, q4, q3, q2 );
    //  cout<<"BC2: Qtot="<<Qtot<<", V="<<V<<endl;
    return V;
}

// a1 -> pi+ pi+ pi- BW
EvtComplex EvtWnPi::BWa( EvtVector4R q )
{
    double const _mA1 = 1.26, _GA1 = 0.4;
    EvtComplex I( 0, 1 );
    double Q2 = q.mass2();
    double GA1 = _GA1 * pi3G( Q2 ) / pi3G( _mA1 * _mA1 );
    EvtComplex denBA1( _mA1 * _mA1 - Q2, -1. * _mA1 * GA1 );
    return _mA1 * _mA1 / denBA1;
}

EvtComplex EvtWnPi::BWf( EvtVector4R q )
{
    double const mf = 0.8, Gf = 0.6;
    EvtComplex I( 0, 1 );
    double Q2 = q.mass2();
    return mf * mf / ( mf * mf - Q2 - I * mf * Gf );
}

EvtComplex EvtWnPi::BWr( EvtVector4R q )
{
    double _mRho = 0.775, _gammaRho = 0.149, _mRhopr = 1.364,
           _gammaRhopr = 0.400, _beta = -0.108;
    double m1 = EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) ),
           m2 = EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) );
    double mQ2 = q.mass2();

    // momenta in the rho->pipi decay
    double dRho = _mRho * _mRho - m1 * m1 - m2 * m2;
    double pPiRho = ( 1.0 / _mRho ) *
                    sqrt( ( dRho * dRho ) / 4.0 - m1 * m1 * m2 * m2 );

    double dRhopr = _mRhopr * _mRhopr - m1 * m1 - m2 * m2;
    double pPiRhopr = ( 1.0 / _mRhopr ) *
                      sqrt( ( dRhopr * dRhopr ) / 4.0 - m1 * m1 * m2 * m2 );

    double dQ = mQ2 - m1 * m1 - m2 * m2;
    double pPiQ = ( 1.0 / sqrt( mQ2 ) ) *
                  sqrt( ( dQ * dQ ) / 4.0 - m1 * m1 * m2 * m2 );

    double gammaRho = _gammaRho * _mRho / sqrt( mQ2 ) *
                      pow( ( pPiQ / pPiRho ), 3 );
    EvtComplex BRhoDem( _mRho * _mRho - mQ2, -1.0 * _mRho * gammaRho );
    EvtComplex BRho = _mRho * _mRho / BRhoDem;

    double gammaRhopr = _gammaRhopr * _mRhopr / sqrt( mQ2 ) *
                        pow( ( pPiQ / pPiRhopr ), 3 );
    EvtComplex BRhoprDem( _mRhopr * _mRhopr - mQ2, -1.0 * _mRho * gammaRhopr );
    EvtComplex BRhopr = _mRhopr * _mRhopr / BRhoprDem;

    return ( BRho + _beta * BRhopr ) / ( 1 + _beta );
}

double EvtWnPi::pi3G( double m2 )
{
    double mPi = EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) );
    double _mRho = 0.775;
    if ( m2 > ( _mRho + mPi ) ) {
        return m2 * ( 1.623 + 10.38 / m2 - 9.32 / ( m2 * m2 ) +
                      0.65 / ( m2 * m2 * m2 ) );
    } else {
        double t1 = m2 - 9.0 * mPi * mPi;
        return 4.1 * pow( t1, 3.0 ) * ( 1.0 - 3.3 * t1 + 5.8 * t1 * t1 );
    };
}

EvtVector4C EvtWnPi::JB( EvtVector4R p1, EvtVector4R p2, EvtVector4R p3,
                         EvtVector4R p4, EvtVector4R p5 )
{
    EvtVector4R Qtot = p1 + p2 + p3 + p4 + p5, Qa = p1 + p2 + p3;
    EvtTensor4C T = ( 1 / Qtot.mass2() ) *
                        EvtGenFunctions::directProd( Qtot, Qtot ) -
                    EvtTensor4C::g();
    EvtVector4R V13 = Qa * ( p2 * ( p1 - p3 ) ) / Qa.mass2() - ( p1 - p3 );
    EvtVector4R V23 = Qa * ( p1 * ( p2 - p3 ) ) / Qa.mass2() - ( p2 - p3 );
    return BWa( Qtot ) * BWa( Qa ) * BWf( p4 + p5 ) *
           ( T.cont1( V13 ) * BWr( p1 + p3 ) + T.cont1( V23 ) * BWr( p2 + p3 ) );
}
