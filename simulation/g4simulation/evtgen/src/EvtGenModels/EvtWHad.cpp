
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

#include "EvtGenModels/EvtWHad.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtTensor4C.hh"

EvtWHad::EvtWHad() : mRho_(), gamma0_(), cK_( 0 )
{
    // cK coefficients from Eur. Phys. J. C39, 41 (2005), arXiv:hep-ph/0409080 [hep-ph]

    // rho(770)
    mRho_.push_back( 0.77511 );
    gamma0_.push_back( 0.1491 );
    cK_.push_back( 1.195 );

    // rho(1450)
    mRho_.push_back( 1.465 );
    gamma0_.push_back( 0.400 );
    cK_.push_back( -0.112 );

    // rho(1700)
    mRho_.push_back( 1.72 );
    gamma0_.push_back( 0.250 );    // rho(1700)
    cK_.push_back( -0.083 );

    // rho(2150), PRD 76 092005
    mRho_.push_back( 2.150 );
    gamma0_.push_back( 0.310 );
    cK_.push_back( 0.0 );
}

EvtComplex EvtWHad::BWKK( double s, int i ) const
{
    double m2 = pow( mRho_[i], 2 );
    EvtComplex qs = pcm( s );
    EvtComplex qm = pcm( m2 );
    if ( abs( qm ) < 1e-10 ) {
        return 0;
    }

    EvtComplex rat = qs / qm;
    EvtComplex rat3 = rat * rat * rat;
    if ( abs( s ) < 1e-10 ) {
        return 0;
    }

    EvtComplex gamma = m2 * rat3 * gamma0_[i] / s;
    EvtComplex I( 0.0, 1.0 );

    EvtComplex denBW = m2 - s - I * sqrt( s ) * gamma;
    if ( abs( denBW ) < 1e-10 ) {
        return 0;
    }

    return cK_[i] * m2 / denBW;
}

EvtVector4C EvtWHad::WCurrent_KSK( const EvtVector4R& pKS,
                                   const EvtVector4R& pKplus ) const
{
    double s = ( pKS + pKplus ).mass2();
    EvtComplex f = BWKK( s, 0 ) + BWKK( s, 1 ) + BWKK( s, 2 );

    return f * ( pKS - pKplus );
}

EvtComplex EvtWHad::pcm( double s ) const
{
    double mpi2 = pow( 0.140, 2 );
    if ( abs( s ) < 1e-10 )
        return 0;

    double pcm2 = 1.0 - 4.0 * mpi2 / s;
    EvtComplex result;

    if ( pcm2 >= 0.0 ) {
        result = EvtComplex( sqrt( pcm2 ), 0.0 );
    } else {
        result = EvtComplex( 0.0, sqrt( -pcm2 ) );
    }

    return result;
}

// =================== W+ -> pi_ current ========================================

EvtVector4C EvtWHad::WCurrent( const EvtVector4R& q1 ) const
{
    return q1;
}

//====================== W+ -> pi+ pi0 current =========================================

EvtVector4C EvtWHad::WCurrent( const EvtVector4R& q1, const EvtVector4R& q2 ) const
{
    return BWr( q1 + q2 ) * ( q1 - q2 );
}

//========================= W+ -> pi+ pi+ pi- current ==============================================

EvtVector4C EvtWHad::WCurrent( const EvtVector4R& q1, const EvtVector4R& q2,
                               const EvtVector4R& q3 ) const
{
    EvtVector4R Q = q1 + q2 + q3;
    EvtVector4R q13 = q1 - q3, q23 = q2 - q3;
    double Q2 = Q.mass2();

    return BWa( Q ) * ( q13 - ( Q * ( Q * q13 ) / Q2 ) * BWr( q2 + q3 ) + q23 -
                        ( Q * ( Q * q23 ) / Q2 ) * BWr( q1 + q3 ) );
}

// ================= W+ -> pi+ pi+ pi- pi- pi+ current with symmetrization ================================

EvtVector4C EvtWHad::WCurrent( const EvtVector4R& q1, const EvtVector4R& q2,
                               const EvtVector4R& q3, const EvtVector4R& q4,
                               const EvtVector4R& q5 ) const
{
    EvtVector4C term1 = JB( q1, q2, q3, q4, q5 );
    EvtVector4C term2 = JB( q5, q2, q3, q4, q1 );
    EvtVector4C term3 = JB( q1, q5, q3, q4, q2 );
    EvtVector4C term4 = JB( q1, q2, q4, q3, q5 );
    EvtVector4C term5 = JB( q5, q2, q4, q3, q1 );
    EvtVector4C term6 = JB( q1, q5, q4, q3, q2 );

    EvtVector4C V = term1 + term2 + term3 + term4 + term5 + term6;
    return V;
}

// =========================W+ -> K+ K- pi+ current ==================================================

EvtVector4C EvtWHad::WCurrent_KKP( const EvtVector4R& pKplus,
                                   const EvtVector4R& pKminus,
                                   const EvtVector4R& pPiPlus ) const
{
    const double mK_892( 0.892 ), gammaK_892( 0.051 );
    const double mA1( 1.239 ), gammaA1( 0.600 );

    EvtVector4R q = pKplus + pKminus + pPiPlus;
    double q2 = q.mass2();
    EvtVector4R pK = pKminus + pPiPlus;
    double pK2 = pK.mass2();

    EvtComplex I( 0.0, 1.0 ), den1, den2;

    den1 = 1.0 / ( q2 - mA1 * mA1 + I * mA1 * gammaA1 );
    den2 = 1.0 / ( pK2 - mK_892 * mK_892 + I * mK_892 * gammaK_892 );

    EvtTensor4C ten = EvtTensor4C::g() -
                      ( 1.0 / q2 ) * EvtGenFunctions::directProd( q, q );

    EvtVector4C vec = den1 * den2 * ( pKminus - pPiPlus );
    vec = ten.cont2( vec );

    return vec;
}

// hadronic current W+ -> K+ pi+ pi-

EvtVector4C EvtWHad::WCurrent_KPP( const EvtVector4R& pKplus,
                                   const EvtVector4R& pPiPlus,
                                   const EvtVector4R& pPiMinus ) const
{
    const double cK1p = 0.210709, cK1r = -0.0152997, cK2p = 0.0945309,
                 cK2r = 0.504315;
    const double mK1_1270 = 1.270, gammaK1_1270 = 0.090, gK1270_Krho = 2.71,
                 gK1270_KsPi = 0.792;
    const double mK1_1400 = 1.400, gammaK1_1400 = 0.174, gK11400_Krho = 0.254,
                 gK11400_KsPi = 2.509;
    const double mK_892 = 0.892, gammaK_892 = 0.051, gK892_Kpi = 3.26;
    const double mRho = 0.770, gammaRho = 0.150, gRho_PiPi = 6.02;

    EvtVector4R q = pKplus + pPiPlus + pPiMinus;
    double q2 = q.mass2(), pp2( 0.0 );

    EvtVector4C curr( 0, 0, 0, 0 ), curr1;
    EvtComplex I( 0.0, 1.0 ), den1, den2;

    // W+ -> K1+(1270) -> K+ rho0 -> K+ pi+ pi-
    den1 = gK1270_Krho /
           ( q2 - mK1_1270 * mK1_1270 + I * mK1_1270 * gammaK1_1270 );
    pp2 = ( pPiPlus + pPiMinus ).mass2();
    den2 = gRho_PiPi / ( pp2 - mRho * mRho + I * mRho * gammaRho );
    curr1 = ( pPiPlus - pPiMinus ) * den1 * den2;
    curr = curr + cK1r * curr1;

    // W+ -> K1+(1270) -> K*(892)0 pi+ -> K+ pi- pi-
    den1 = gK1270_KsPi /
           ( q2 - mK1_1270 * mK1_1270 + I * mK1_1270 * gammaK1_1270 );
    pp2 = ( pKplus + pPiMinus ).mass2();
    den2 = gK892_Kpi / ( pp2 - mK_892 * mK_892 + I * mK_892 * gammaK_892 );
    curr1 = ( pKplus - pPiMinus ) * den1 * den2;
    curr = curr + cK1p * curr1;

    // W+ -> K1+(1400) -> K+ rho0 -> K+ pi+ pi-
    den1 = gK11400_Krho /
           ( q2 - mK1_1400 * mK1_1400 + I * mK1_1400 * gammaK1_1400 );
    pp2 = ( pPiMinus + pPiPlus ).mass2();
    den2 = gRho_PiPi / ( pp2 - mRho * mRho + I * mRho * gammaRho );
    curr1 = ( pPiPlus - pPiMinus ) * den1 * den2;
    curr = curr + cK2r * curr1;

    // W+ -> K1+(1400) -> K*(892)0 pi+ -> K+ pi- pi+
    den1 = gK11400_KsPi /
           ( q2 - mK1_1400 * mK1_1400 + I * mK1_1400 * gammaK1_1400 );
    pp2 = ( pKplus + pPiMinus ).mass2();
    den2 = gK892_Kpi / ( pp2 - mK_892 * mK_892 + I * mK_892 * gammaK_892 );
    curr1 = ( pKplus - pPiPlus ) * den1 * den2;
    curr = curr + cK2p * curr1;

    EvtTensor4C ten = EvtTensor4C::g() -
                      ( 1.0 / q2 ) * EvtGenFunctions::directProd( q, q );
    curr = ten.cont2( curr );

    return curr;
}

// a1 -> pi+ pi+ pi- BW

EvtComplex EvtWHad::BWa( const EvtVector4R& q ) const
{
    double _mA1( 1.26 ), _GA1( 0.4 );
    double _mA1Sq = _mA1 * _mA1;
    EvtComplex I( 0.0, 1.0 );
    double Q2 = q.mass2();
    double GA1 = _GA1 * pi3G( Q2 ) / pi3G( _mA1Sq );

    EvtComplex denBA1( _mA1Sq - Q2, -1.0 * _mA1 * GA1 );

    return _mA1Sq / denBA1;
}

EvtComplex EvtWHad::BWf( const EvtVector4R& q ) const
{
    double mf( 0.8 ), Gf( 0.6 );
    double mfSq = mf * mf;
    EvtComplex I( 0.0, 1.0 );
    double Q2 = q.mass2();
    return mfSq / ( mfSq - Q2 - I * mf * Gf );
}

EvtComplex EvtWHad::BWr( const EvtVector4R& q ) const
{
    double _mRho( 0.775 ), _gammaRho( 0.149 ), _mRhopr( 1.364 ),
        _gammaRhopr( 0.400 ), _beta( -0.108 );
    double m1 = EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) ),
           m2 = EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) );
    double mQ2 = q.mass2();

    double m1Sq = m1 * m1, m2Sq = m2 * m2, m12Sq = m1Sq * m2Sq;

    // momenta in the rho->pipi decay
    double dRho = _mRho * _mRho - m1Sq - m2Sq;
    double pPiRho = ( 1.0 / _mRho ) * sqrt( ( dRho * dRho ) / 4.0 - m12Sq );

    double dRhopr = _mRhopr * _mRhopr - m1Sq - m2Sq;
    double pPiRhopr = ( 1.0 / _mRhopr ) *
                      sqrt( ( dRhopr * dRhopr ) / 4.0 - m12Sq );

    double dQ = mQ2 - m1Sq - m2Sq;
    double pPiQ = ( 1.0 / sqrt( mQ2 ) ) * sqrt( ( dQ * dQ ) / 4.0 - m12Sq );

    double gammaRho = _gammaRho * _mRho / sqrt( mQ2 ) *
                      pow( ( pPiQ / pPiRho ), 3 );
    EvtComplex BRhoDem( _mRho * _mRho - mQ2, -1.0 * _mRho * gammaRho );
    EvtComplex BRho = _mRho * _mRho / BRhoDem;

    double gammaRhopr = _gammaRhopr * _mRhopr / sqrt( mQ2 ) *
                        pow( ( pPiQ / pPiRhopr ), 3 );
    EvtComplex BRhoprDem( _mRhopr * _mRhopr - mQ2, -1.0 * _mRho * gammaRhopr );
    EvtComplex BRhopr = _mRhopr * _mRhopr / BRhoprDem;

    return ( BRho + _beta * BRhopr ) / ( 1.0 + _beta );
}

double EvtWHad::pi3G( double m2 ) const
{
    double mPi = EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) );
    double mRho( 0.775 );
    double val( 0.0 );

    if ( m2 > ( mRho + mPi ) ) {
        val = m2 * ( 1.623 + 10.38 / m2 - 9.32 / ( m2 * m2 ) +
                     0.65 / ( m2 * m2 * m2 ) );
    } else {
        double t1 = m2 - 9.0 * mPi * mPi;
        val = 4.1 * pow( t1, 3.0 ) * ( 1.0 - 3.3 * t1 + 5.8 * t1 * t1 );
    }

    return val;
}

EvtVector4C EvtWHad::JB( const EvtVector4R& p1, const EvtVector4R& p2,
                         const EvtVector4R& p3, const EvtVector4R& p4,
                         const EvtVector4R& p5 ) const
{
    EvtVector4R Qtot = p1 + p2 + p3 + p4 + p5, Qa = p1 + p2 + p3;
    EvtTensor4C T = ( 1.0 / Qtot.mass2() ) *
                        EvtGenFunctions::directProd( Qtot, Qtot ) -
                    EvtTensor4C::g();

    EvtVector4R p13 = p1 - p3, p23 = p2 - p3;
    EvtVector4R V13 = Qa * ( p2 * p13 ) / Qa.mass2() - p13;
    EvtVector4R V23 = Qa * ( p1 * p23 ) / Qa.mass2() - p23;

    return BWa( Qtot ) * BWa( Qa ) * BWf( p4 + p5 ) *
           ( T.cont1( V13 ) * BWr( p1 + p3 ) + T.cont1( V23 ) * BWr( p2 + p3 ) );
}
