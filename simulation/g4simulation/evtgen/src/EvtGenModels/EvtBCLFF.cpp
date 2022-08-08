
/***********************************************************************
* Copyright 1998-2021 CERN for the benefit of the EvtGen authors       *
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

#include "EvtGenModels/EvtBCLFF.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <cmath>
#include <cstdlib>
#include <string>

EvtBCLFF::EvtBCLFF( int numarg, double* arglist ) :
    m_numBCLFFCoefficients( numarg )
{
    if ( numarg > 19 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Too large number of coefficients!\n";
        ::abort();
    }
    for ( int i = 0; i < m_numBCLFFCoefficients; ++i ) {
        m_BCLFFCoefficients[i] = arglist[i];
    }
}

void EvtBCLFF::getscalarff( EvtId parent, EvtId daughter, double t, double,
                            double* fpf, double* f0f )
{
    if ( m_numBCLFFCoefficients != 8 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Wrong number of arguments for EvtBCLFF::getscalarff!\n";
        ::abort();
    }

    const auto mB = EvtPDL::getMeanMass( parent );
    const auto mM = EvtPDL::getMeanMass( daughter );

    const auto tplus = ( mB + mM ) * ( mB + mM );
    const auto tzero = ( mB + mM ) * ( std::sqrt( mB ) - std::sqrt( mM ) ) *
                       ( std::sqrt( mB ) - std::sqrt( mM ) );

    const auto mR2 = m_resonance1Minus * m_resonance1Minus;
    const auto pole = 1.0 / ( 1.0 - t / mR2 );

    const std::array<double, 4> bplus{ m_BCLFFCoefficients[0],
                                       m_BCLFFCoefficients[1],
                                       m_BCLFFCoefficients[2],
                                       m_BCLFFCoefficients[3] };
    const std::array<double, 4> bzero{ m_BCLFFCoefficients[4],
                                       m_BCLFFCoefficients[5],
                                       m_BCLFFCoefficients[6],
                                       m_BCLFFCoefficients[7] };

    const auto N_fpf = bplus.size();
    const auto N_f0f = bzero.size();

    auto z = [tplus, tzero]( decltype( t ) q2 ) {
        const auto term1 = std::sqrt( tplus - q2 );
        const auto term2 = std::sqrt( tplus - tzero );
        return ( term1 - term2 ) / ( term1 + term2 );
    };

    double sum_fpf = 0;
    for ( unsigned int n = 0; n < N_fpf; ++n ) {
        sum_fpf += bplus[n] * ( std::pow( z( t ), n ) -
                                std::pow( -1, n - N_fpf ) * n / N_fpf *
                                    std::pow( z( t ), N_fpf ) );
    }
    *fpf = pole * sum_fpf;

    double sum_f0f = 0;
    for ( unsigned int n = 0; n < N_f0f; ++n ) {
        sum_f0f += bzero[n] * std::pow( z( t ), n );
    }
    *f0f = sum_f0f;
}

void EvtBCLFF::getvectorff( EvtId parent, EvtId daughter, double t, double,
                            double* a1f, double* a2f, double* vf, double* a0f )
{
    if ( m_numBCLFFCoefficients != 11 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Wrong number of arguments for EvtBCLFF::getvectorff!\n";
        ::abort();
    }

    const auto mB = EvtPDL::getMeanMass( parent );
    const auto mB2 = mB * mB;
    const auto mM = EvtPDL::getMeanMass( daughter );
    const auto mM2 = mM * mM;

    const auto tplus = ( mB + mM ) * ( mB + mM );
    const auto tminus = ( mB - mM ) * ( mB - mM );
    const auto tzero = tplus * ( 1.0 - std::sqrt( 1.0 - tminus / tplus ) );

    const auto mR2A0 = m_resonance0Minus * m_resonance0Minus;
    const auto mR2A1 = m_resonance1Plus * m_resonance1Plus;
    const auto mR2A12 = m_resonance1Plus * m_resonance1Plus;
    const auto mR2V = m_resonance1Minus * m_resonance1Minus;

    const auto poleA0 = 1.0 / ( 1.0 - t / mR2A0 );
    const auto poleA1 = 1.0 / ( 1.0 - t / mR2A1 );
    const auto poleA12 = 1.0 / ( 1.0 - t / mR2A12 );
    const auto poleV = 1.0 / ( 1.0 - t / mR2V );

    const std::array<double, 3> A0{
        8 * mB * mM / ( mB2 - mM2 ) * m_BCLFFCoefficients[5],
        m_BCLFFCoefficients[0], m_BCLFFCoefficients[1] };
    const std::array<double, 3> A1{ m_BCLFFCoefficients[2],
                                    m_BCLFFCoefficients[3],
                                    m_BCLFFCoefficients[4] };
    const std::array<double, 3> A12{ m_BCLFFCoefficients[5],
                                     m_BCLFFCoefficients[6],
                                     m_BCLFFCoefficients[7] };
    const std::array<double, 3> V{ m_BCLFFCoefficients[8], m_BCLFFCoefficients[9],
                                   m_BCLFFCoefficients[10] };

    auto z = [tplus, tzero]( decltype( t ) q2 ) {
        const auto term1 = std::sqrt( tplus - q2 );
        const auto term2 = std::sqrt( tplus - tzero );
        return ( term1 - term2 ) / ( term1 + term2 );
    };

    auto sum = [&z]( decltype( t ) q2, std::array<double, 3> par ) {
        double tot = 0.0;
        for ( unsigned int n = 0; n < par.size(); ++n ) {
            tot += par[n] * std::pow( z( q2 ) - z( 0.0 ), n );
        }
        return tot;
    };

    auto kaellen = [mB, mM]( decltype( t ) q2 ) {
        return ( ( mB + mM ) * ( mB + mM ) - q2 ) *
               ( ( mB - mM ) * ( mB - mM ) - q2 );
    };

    const auto ffA0 = poleA0 * sum( t, A0 );
    const auto ffA1 = poleA1 * sum( t, A1 );
    const auto ffA12 = poleA12 * sum( t, A12 );
    const auto ffV = poleV * sum( t, V );

    const auto ffA2 = ( ( mB + mM ) * ( mB + mM ) * ( mB2 - mM2 - t ) * ffA1 -
                        ( 16 * mB * mM2 * ( mB + mM ) ) * ffA12 ) /
                      kaellen( t );

    *a0f = ffA0;
    *a1f = ffA1;
    *a2f = ffA2;
    *vf = ffV;
}

void EvtBCLFF::gettensorff( EvtId, EvtId, double, double, double*, double*,
                            double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :gettensorff in EvtBCLFF.\n";
    ::abort();
}

void EvtBCLFF::getbaryonff( EvtId, EvtId, double, double, double*, double*,
                            double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtBCLFF.\n";
    ::abort();
}

void EvtBCLFF::getdiracff( EvtId, EvtId, double, double, double*, double*,
                           double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtBCLFF.\n";
    ::abort();
}

void EvtBCLFF::getraritaff( EvtId, EvtId, double, double, double*, double*,
                            double*, double*, double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtBCLFF.\n";
    ::abort();
}
