
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

#include "EvtGenModels/EvtSVPHelAmp.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor3C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector3C.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtSVPHelAmp::getName()
{
    return "SVP_HELAMP";
}

EvtDecayBase* EvtSVPHelAmp::clone()
{
    return new EvtSVPHelAmp;
}

void EvtSVPHelAmp::initProbMax()
{
    setProbMax( 2.0 * ( getArg( 0 ) * getArg( 0 ) + getArg( 2 ) * getArg( 2 ) ) );
}

void EvtSVPHelAmp::init()
{
    // check that there are 4 arguments
    checkNArg( 4 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::PHOTON );
}

void EvtSVPHelAmp::decay( EvtParticle* p )
{
    SVPHel( p, _amp2, getDaug( 0 ), getDaug( 1 ),
            EvtComplex( getArg( 0 ) * cos( getArg( 1 ) ),
                        getArg( 0 ) * sin( getArg( 1 ) ) ),
            EvtComplex( getArg( 2 ) * cos( getArg( 3 ) ),
                        getArg( 2 ) * sin( getArg( 3 ) ) ) );

    return;
}

void EvtSVPHelAmp::SVPHel( EvtParticle* parent, EvtAmp& amp, EvtId n_v1,
                           EvtId n_v2, const EvtComplex& hp, const EvtComplex& hm )
{
    //  Routine to decay a vector into a vector and scalar.  Started
    //  by ryd on Oct 17, 1996.

    // This routine is adopted from EvtSVVHel and since there is
    // a photon that can not have helicity 0 this is put in by
    // setting the h0 amplitude to 0.
    EvtComplex h0 = EvtComplex( 0.0, 0.0 );

    int tndaug = 2;
    EvtId tdaug[2];
    tdaug[0] = n_v1;
    tdaug[1] = n_v2;

    parent->initializePhaseSpace( tndaug, tdaug );

    EvtParticle *v1, *v2;
    v1 = parent->getDaug( 0 );
    v2 = parent->getDaug( 1 );

    EvtVector4R momv1 = v1->getP4();

    EvtVector3R v1dir( momv1.get( 1 ), momv1.get( 2 ), momv1.get( 3 ) );
    v1dir = v1dir / v1dir.d3mag();

    EvtComplex a = -0.5 * ( hp + hm );
    EvtComplex b = EvtComplex( 0.0, 0.5 ) * ( hp - hm );
    EvtComplex c = h0 + 0.5 * ( hp + hm );

    EvtTensor3C M = a * EvtTensor3C::id() + b * EvtGenFunctions::eps( v1dir ) +
                    c * EvtGenFunctions::directProd( v1dir, v1dir );

    EvtVector3C t0 = M.cont1( v1->eps( 0 ).vec().conj() );
    EvtVector3C t1 = M.cont1( v1->eps( 1 ).vec().conj() );
    EvtVector3C t2 = M.cont1( v1->eps( 2 ).vec().conj() );

    EvtVector3C eps0 = v2->epsParentPhoton( 0 ).vec().conj();
    EvtVector3C eps1 = v2->epsParentPhoton( 1 ).vec().conj();

    amp.vertex( 0, 0, t0 * eps0 );
    amp.vertex( 0, 1, t0 * eps1 );
    //amp.vertex(0,2,t1*eps0*0.);

    amp.vertex( 1, 0, t1 * eps0 );
    amp.vertex( 1, 1, t1 * eps1 );
    //amp.vertex(1,2,t1*eps0*0.);

    amp.vertex( 2, 0, t2 * eps0 );
    amp.vertex( 2, 1, t2 * eps1 );
    //amp.vertex(2,2,t1*eps0*0.);

    return;
}
