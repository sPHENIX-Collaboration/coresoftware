
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

#include "EvtGenBase/EvtPhotonParticle.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <iostream>
#include <math.h>
#include <stdlib.h>
using std::endl;

void EvtPhotonParticle::init( EvtId part_n, const EvtVector4R& p4 )
{
    init( part_n, p4.get( 0 ), p4.get( 1 ), p4.get( 2 ), p4.get( 3 ) );
}

void EvtPhotonParticle::init( EvtId part_n, double e, double px, double py,
                              double pz )
{
    _validP4 = true;
    setp( e, px, py, pz );
    setpart_num( part_n );

    setLifetime();

    //defere calculation of basis vectors untill they are needed!
    _evalBasis = 0;
}

EvtVector4C EvtPhotonParticle::epsParentPhoton( int i )
{
    if ( !_evalBasis ) {
        _evalBasis = 1;
        eps1.set( EvtComplex( 0.0, 0.0 ), EvtComplex( -1.0 / sqrt( 2.0 ), 0.0 ),
                  EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), EvtComplex( 0.0, 0.0 ) );
        eps2.set( EvtComplex( 0.0, 0.0 ), EvtComplex( 1.0 / sqrt( 2.0 ), 0.0 ),
                  EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), EvtComplex( 0.0, 0.0 ) );

        // These are for photon along z axis.  Rotate to get
        // correct direction...

        double phi, theta;

        EvtVector4R p = this->getP4();

        double px = p.get( 1 );
        double py = p.get( 2 );
        double pz = p.get( 3 );

        phi = atan2( py, px );
        theta = acos( pz / sqrt( px * px + py * py + pz * pz ) );
        eps1.applyRotateEuler( phi, theta, -phi );
        eps2.applyRotateEuler( phi, theta, -phi );
    }

    EvtVector4C temp;

    switch ( i ) {
        case 0:
            temp = eps1;
            break;
        case 1:
            temp = eps2;
            break;
        default:
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "EvtPhotonParticle.cc: Asked "
                << "for state:" << i << endl;
            ::abort();
            break;
    }

    return temp;
}

EvtVector4C EvtPhotonParticle::epsPhoton( int )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "EvtPhotonParticle.cc: Can not get "
        << "state in photons restframe." << endl;
    ;
    ::abort();
    return EvtVector4C();
}

EvtSpinDensity EvtPhotonParticle::rotateToHelicityBasis() const
{
    EvtVector4C eplus( 0.0, -1.0 / sqrt( 2.0 ),
                       EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), 0.0 );
    EvtVector4C eminus( 0.0, 1.0 / sqrt( 2.0 ),
                        EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), 0.0 );

    //Really uggly have to cast away constness because the
    //function epsParentPhoton caches the state vectors...
    EvtVector4C e1 = ( (EvtParticle*)this )->epsParentPhoton( 0 );
    EvtVector4C e2 = ( (EvtParticle*)this )->epsParentPhoton( 1 );

    EvtSpinDensity R;
    R.setDim( 2 );

    R.set( 0, 0, ( eplus.conj() ) * e1 );
    R.set( 0, 1, ( eplus.conj() ) * e2 );

    R.set( 1, 0, ( eminus.conj() ) * e1 );
    R.set( 1, 1, ( eminus.conj() ) * e2 );

    return R;
}

EvtSpinDensity EvtPhotonParticle::rotateToHelicityBasis( double alpha,
                                                         double beta,
                                                         double gamma ) const
{
    EvtVector4C eplus( 0.0, -1.0 / sqrt( 2.0 ),
                       EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), 0.0 );
    EvtVector4C eminus( 0.0, 1.0 / sqrt( 2.0 ),
                        EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), 0.0 );

    eplus.applyRotateEuler( alpha, beta, gamma );
    eminus.applyRotateEuler( alpha, beta, gamma );

    //Really uggly have to cast away constness because the
    //function epsParentPhoton caches the state vectors...
    EvtVector4C e1 = ( (EvtParticle*)this )->epsParentPhoton( 0 );
    EvtVector4C e2 = ( (EvtParticle*)this )->epsParentPhoton( 1 );

    EvtSpinDensity R;
    R.setDim( 2 );

    R.set( 0, 0, ( eplus.conj() ) * e1 );
    R.set( 0, 1, ( eplus.conj() ) * e2 );

    R.set( 1, 0, ( eminus.conj() ) * e1 );
    R.set( 1, 1, ( eminus.conj() ) * e2 );

    return R;
}
