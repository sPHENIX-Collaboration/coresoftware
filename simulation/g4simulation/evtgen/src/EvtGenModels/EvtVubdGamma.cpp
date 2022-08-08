
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

#include "EvtGenModels/EvtVubdGamma.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDiLog.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <math.h>

//----------------
// Constructors --
//----------------

EvtVubdGamma::EvtVubdGamma( const double& alphas )
{
    _alphas = alphas;

    // the range for the delta distribution in p2 is from _epsilon1 to
    // _epsilon2. It was checked with the single differential formulae
    // in the paper that these values are small enough to imitate p2 = 0
    // for the regular terms.
    // The ()* distributions, however need further treatment. In order to
    // generate the correct spectrum in z a threshold need to be computed
    // from the desired value of the coupling alphas. The idea is that
    // for z=1 p2=0 is not allowed and therefore the part of dGamma proportional
    // to delta(p2) should go to 0 for z->1.
    // Using equation (3.1) and (3.2) it is possible to find the correct value
    // for log(_epsilon3) from this requirement.

    _epsilon1 = 1e-10;
    _epsilon2 = 1e-5;
    if ( alphas > 0 ) {
        double lne3 = 9. / 16. - 2 * EvtConst::pi * EvtConst::pi / 3. +
                      6 * EvtConst::pi / 4 / alphas;
        if ( lne3 > 0 )
            lne3 = -7. / 4. - sqrt( lne3 );
        else
            lne3 = -7. / 4.;
        _epsilon3 = exp( lne3 );
    } else
        _epsilon3 = 1;
}

//-----------
// Methods --
//-----------

double EvtVubdGamma::getdGdxdzdp( const double& x, const double& z,
                                  const double& p2 )
{
    // check phase space

    double xb = ( 1 - x );

    if ( x < 0 || x > 1 || z < xb || z > ( 1 + xb ) )
        return 0;

    double p2min = ( 0 > z - 1. ? 0 : z - 1. );
    double p2max = ( 1. - x ) * ( z - 1. + x );

    if ( p2 < p2min || p2 > p2max )
        return 0;

    //  // check the phase space
    //  return 1.;

    double dG;

    if ( p2 > _epsilon1 && p2 < _epsilon2 ) {
        double W1 = getW1delta( x, z );
        double W4plus5 = getW4plus5delta( x, z );

        dG = 12. * delta( p2, p2min, p2max ) *
             ( ( 1. + xb - z ) * ( z - xb ) * W1 + xb * ( z - xb ) * ( W4plus5 ) );
    } else {
        double W1 = getW1nodelta( x, z, p2 );
        double W2 = getW2nodelta( x, z, p2 );
        double W3 = getW3nodelta( x, z, p2 );
        double W4 = getW4nodelta( x, z, p2 );
        double W5 = getW5nodelta( x, z, p2 );

        dG = 12. *
             ( ( 1. + xb - z ) * ( z - xb - p2 ) * W1 + ( 1. - z + p2 ) * W2 +
               ( xb * ( z - xb ) - p2 ) * ( W3 + W4 + W5 ) );
    }
    return dG;
}

double EvtVubdGamma::delta( const double& x, const double& xmin,
                            const double& xmax )
{
    if ( xmin > 0 || xmax < 0 )
        return 0.;
    if ( _epsilon1 < x && x < _epsilon2 )
        return 1. / ( _epsilon2 - _epsilon1 );
    return 0.0;
}

double EvtVubdGamma::getW1delta( const double&, const double& z )
{
    double mz = 1. - z;

    double lz;
    if ( z == 1 )
        lz = -1.;
    else
        lz = log( z ) / ( 1. - z );

    // ddilog_(&z) is actually the dilog of (1-z) in maple,
    // also in Neuberts paper the limit dilog(1) = pi^2/6 is used
    // this corresponds to maple's dilog(0), so
    // I take ddilog_(&mz) where mz=1-z in order to satisfy Neubert's definition
    // and to compare with Maple the argument in maple should be (1-mz) ...

    double dl = 4. * EvtDiLog::DiLog( mz ) + 4. * pow( EvtConst::pi, 2 ) / 3.;

    double w = -( 8. * pow( log( z ), 2 ) - 10. * log( z ) + 2. * lz + dl + 5. ) +
               ( 8. * log( z ) - 7. ) * log( _epsilon3 ) -
               2. * pow( log( _epsilon3 ), 2 );

    return ( 1. + w * _alphas / 3. / EvtConst::pi );
}

double EvtVubdGamma::getW1nodelta( const double&, const double& z,
                                   const double& p2 )
{
    double z2 = z * z;
    double t2 = 1. - 4. * p2 / z2;
    double t = sqrt( t2 );

    double w = 0;
    if ( p2 > _epsilon2 )
        w += 4. / p2 * ( log( ( 1. + t ) / ( 1. - t ) ) / t + log( p2 / z2 ) ) +
             1. - ( 8. - z ) * ( 2. - z ) / z2 / t2 +
             ( ( 2. - z ) / 2. / z + ( 8. - z ) * ( 2. - z ) / 2. / z2 / t2 ) *
                 log( ( 1. + t ) / ( 1. - t ) ) / t;
    if ( p2 > _epsilon3 )
        w += ( 8. * log( z ) - 7. ) / p2 - 4. * log( p2 ) / p2;

    return w * _alphas / 3. / EvtConst::pi;
}

double EvtVubdGamma::getW2nodelta( const double&, const double& z,
                                   const double& p2 )
{
    double z2 = z * z;
    double t2 = 1. - 4. * p2 / z2;
    double t = sqrt( t2 );
    double w11 = ( 32. - 8. * z + z2 ) / 4. / z / t2;

    double w = 0;
    if ( p2 > _epsilon2 )
        w -= ( z * t2 / 8. + ( 4. - z ) / 4. + w11 / 2. ) *
             log( ( 1. + t ) / ( 1. - t ) ) / t;
    if ( p2 > _epsilon2 )
        w += ( 8. - z ) / 4. + w11;

    return ( w * _alphas / 3. / EvtConst::pi );
}

double EvtVubdGamma::getW3nodelta( const double&, const double& z,
                                   const double& p2 )
{
    double z2 = z * z;
    double t2 = 1. - 4. * p2 / z2;
    double t4 = t2 * t2;
    double t = sqrt( t2 );

    double w = 0;

    if ( p2 > _epsilon2 )
        w += ( z * t2 / 16. + 5. * ( 4. - z ) / 16. -
               ( 64. + 56. * z - 7. * z2 ) / 16. / z / t2 +
               3. * ( 12. - z ) / 16. / t4 ) *
             log( ( 1. + t ) / ( 1. - t ) ) / t;
    if ( p2 > _epsilon2 )
        w += -( 8. - 3. * z ) / 8. + ( 32. + 22. * z - 3. * z2 ) / 4. / z / t2 -
             3. * ( 12. - z ) / 8. / t4;

    return ( w * _alphas / 3. / EvtConst::pi );
}

double EvtVubdGamma::getW4nodelta( const double&, const double& z,
                                   const double& p2 )
{
    double z2 = z * z;
    double t2 = 1. - 4. * p2 / z2;
    double t4 = t2 * t2;
    double t = sqrt( t2 );

    double w = 0;

    if ( p2 > _epsilon2 )
        w -= ( ( 8. - 3. * z ) / 4. / z - ( 22. - 3. * z ) / 2. / z / t2 +
               3. * ( 12. - z ) / 4. / z / t4 ) *
             log( ( 1. + t ) / ( 1. - t ) ) / t;
    if ( p2 > _epsilon2 )
        w += -1. - ( 32. - 5. * z ) / 2. / z / t2 +
             3. * ( 12. - z ) / 2. / z / t4;

    return w * _alphas / 3. / EvtConst::pi;
}

double EvtVubdGamma::getW4plus5delta( const double&, const double& z )
{
    double w = 0;

    if ( z == 1 )
        w = -2;
    else
        w = 2. * log( z ) / ( 1. - z );

    return ( w * _alphas / 3. / EvtConst::pi );
}

double EvtVubdGamma::getW5nodelta( const double&, const double& z,
                                   const double& p2 )
{
    double z2 = z * z;
    double t2 = 1. - 4. * p2 / z2;
    double t4 = t2 * t2;
    double t = sqrt( t2 );

    double w = 0;
    if ( p2 > _epsilon2 )
        w += ( 1. / 4. / z - ( 2. - z ) / 2. / z2 / t2 +
               3. * ( 12. - z ) / 4. / z2 / t4 ) *
             log( ( 1. + t ) / ( 1. - t ) ) / t;
    if ( p2 > _epsilon2 )
        w += -( 8. + z ) / 2. / z2 / t2 - 3. * ( 12. - z ) / 2. / z2 / t4;

    return ( w * _alphas / 3. / EvtConst::pi );
}
