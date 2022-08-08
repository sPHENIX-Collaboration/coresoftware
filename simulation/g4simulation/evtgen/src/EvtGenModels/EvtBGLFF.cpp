
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

#include "EvtGenModels/EvtBGLFF.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"

#include <cmath>
#include <cstdlib>

// BGL (N=3) for scalar meson i.e. B->Dlv  (l=e,mu,tau)
EvtBGLFF::EvtBGLFF( double bglap_0, double bglap_1, double bglap_2,
                    double bglap_3, double bgla0_0, double bgla0_1,
                    double bgla0_2, double bgla0_3 ) :
    m_ap_0( bglap_0 ),
    m_ap_1( bglap_1 ),
    m_ap_2( bglap_2 ),
    m_ap_3( bglap_3 ),
    m_a0_0( bgla0_0 ),
    m_a0_1( bgla0_1 ),
    m_a0_2( bgla0_2 ),
    m_a0_3( bgla0_3 )
{
}

// BGL for vector meson i.e. B->D*lv (l=e,mu), and should not be used to taus.
EvtBGLFF::EvtBGLFF( double bgla_0, double bgla_1, double bglb_0, double bglb_1,
                    double bglc_1, double bglc_2 ) :
    m_a_0( bgla_0 ),
    m_a_1( bgla_1 ),
    m_b_0( bglb_0 ),
    m_b_1( bglb_1 ),
    m_c_1( bglc_1 ),
    m_c_2( bglc_2 )
{
}

// Use dispersion relation parametrization from
// C.G.Boyd, B.Grinstein, R.F.Lebed, Phys. Rev. Lett. 74,4603(1995)
// and
// R.Glattauer, etc. (Belle) Phys. Rev. D 93,032006 (2016)
void EvtBGLFF::getscalarff( EvtId parent, EvtId, double t, double mass,
                            double* fp, double* f0 )
{
    const double mb = EvtPDL::getMeanMass( parent );
    const double r = mass / mb;
    const double w = ( ( mb * mb ) + ( mass * mass ) - t ) / ( 2.0 * mb * mass );
    const double z = ( std::sqrt( w + 1.0 ) - sqrt( 2. ) ) /
                     ( std::sqrt( w + 1.0 ) + std::sqrt( 2. ) );
    const double p_i = 1.0;

    const double phi_sub = ( 1.0 + r ) * ( 1.0 - z ) +
                           2.0 * std::sqrt( r ) * ( 1.0 + z );
    const double g_sub = ( 4.0 * r ) / std::pow( 1.0 + r, 2 );
    const double phi_p = 1.1213 * std::pow( 1.0 + z, 2. ) * sqrt( 1.0 - z ) *
                         std::pow( phi_sub, -5 );
    const double phi_0 = 0.5299 * ( 1.0 + z ) * std::pow( 1.0 - z, 1.5 ) *
                         std::pow( phi_sub, -4 );

    *fp = g_sub * ( m_ap_0 + m_ap_1 * z + m_ap_2 * z * z + m_ap_3 * z * z * z ) /
          ( p_i * phi_p );
    *f0 = g_sub * ( m_a0_0 + m_a0_1 * z + m_a0_2 * z * z + m_a0_3 * z * z * z ) /
          ( p_i * phi_0 );

    return;
}

// Use z expansion parametrization from
// C.G.Boyd, B.Grinstein and R.F.Lebed, Phys. Rev. D 56,6895(1997)
// and
// B.Grinstein, A.Kobach, Phys. Lett. B 771(2017)359-364
void EvtBGLFF::getvectorff( EvtId parent, EvtId, double t, double mass,
                            double* a1f, double* a2f, double* vf, double* a0f )
{
    double mb = EvtPDL::getMeanMass( parent );
    double w = ( ( mb * mb ) + ( mass * mass ) - t ) / ( 2. * mb * mass );

    // Form factors have a general form, with parameters passed in
    // from the arguments.

    const double r = mass / mb;
    const double z = ( std::sqrt( w + 1. ) - std::sqrt( 2. ) ) /
                     ( std::sqrt( w + 1 ) + std::sqrt( 2. ) );
    const double rstar = ( 2. * std::sqrt( mb * mass ) ) / ( mb + mass );
    constexpr double chiT_plus33 = 5.28e-4;
    constexpr double chiT_minus33 = 3.07e-4;
    constexpr double n_i = 2.6;
    constexpr double axialvector_poles[4]{ 6.730, 6.736, 7.135, 7.142 };
    constexpr double vector_poles[4]{ 6.337, 6.899, 7.012, 7.280 };

    const double c_0 = ( mb - mass ) / mb * std::sqrt( 0.5 ) /
                       ( 1.0 + r + 2. * std::sqrt( r ) ) * m_b_0;

    const double phi_g = std::sqrt( 256. * n_i / ( 3. * M_PI * chiT_plus33 ) ) *
                         r * r * ( 1. + z ) * ( 1. + z ) / std::sqrt( 1. - z ) /
                         std::pow( ( 1. + r ) * ( 1. - z ) +
                                       2. * std::sqrt( r ) * ( 1. + z ),
                                   4. );
    const double phi_f =
        1. / ( mb * mb ) * std::sqrt( 16. * n_i / ( 3. * M_PI * chiT_minus33 ) ) *
        r * ( 1. + z ) * std::pow( 1. - z, 1.5 ) /
        std::pow( ( 1. + r ) * ( 1. - z ) + 2. * std::sqrt( r ) * ( 1. + z ), 4. );
    const double phi_F1 =
        1. / ( mb * mb * mb ) *
        std::sqrt( 8. * n_i / ( 3. * M_PI * chiT_minus33 ) ) * r * ( 1. + z ) *
        std::pow( 1. - z, 2.5 ) /
        std::pow( ( 1. + r ) * ( 1. - z ) + 2 * std::sqrt( r ) * ( 1. + z ), 5. );

    double p_g = 1.;
    double p_f = 1.;
    const double term3 = std::sqrt( ( mb + mass ) * ( mb + mass ) -
                                    ( mb - mass ) * ( mb - mass ) );

    for ( int i = 0; i < 4; i++ ) {
        const double term1 = std::sqrt( ( mb + mass ) * ( mb + mass ) -
                                        vector_poles[i] * vector_poles[i] );
        const double term2 = std::sqrt( ( mb + mass ) * ( mb + mass ) -
                                        axialvector_poles[i] *
                                            axialvector_poles[i] );
        const double z_p1 = ( term1 - term3 ) / ( term1 + term3 );
        p_g = p_g * ( z - z_p1 ) / ( 1 - z * z_p1 );
        const double z_p2 = ( term2 - term3 ) / ( term2 + term3 );
        p_f = p_f * ( z - z_p2 ) / ( 1 - z * z_p2 );
    }

    const double g = 1. / p_g / phi_g * ( m_a_0 + m_a_1 * z );
    const double f = 1. / p_f / phi_f * ( m_b_0 + m_b_1 * z );
    const double F1 = 1. / p_f / phi_F1 * ( c_0 + m_c_1 * z + m_c_2 * z * z );

    const double ha1 = f / std::sqrt( mb * mass ) / ( 1. + w );
    const double r1 = ( w + 1. ) * mb * mass * g / f;
    const double r2 = ( w - r ) / ( w - 1.0 ) - F1 / mb / ( w - 1.0 ) / f;

    *a1f = ( w + 1. ) / 2. * rstar * ha1;
    *a2f = ( r2 / rstar ) * ha1;
    *vf = ( r1 / rstar ) * ha1;
    *a0f = 0;    // a0f is related to B->D* tau nu decay.
    // The class should not be used for taus, due to the lack of fitted parameters in a0f amplitude.

    return;
}

void EvtBGLFF::gettensorff( EvtId, EvtId, double, double, double*, double*,
                            double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :gettensorff in EvtBGLFF.\n";
    ::abort();
}

void EvtBGLFF::getbaryonff( EvtId, EvtId, double, double, double*, double*,
                            double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbayronff in EvtBGLFF.\n";
    ::abort();
}

void EvtBGLFF::getdiracff( EvtId, EvtId, double, double, double*, double*,
                           double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtBGLFF.\n";
    ::abort();
}

void EvtBGLFF::getraritaff( EvtId, EvtId, double, double, double*, double*,
                            double*, double*, double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtBGLFF.\n";
    ::abort();
}
