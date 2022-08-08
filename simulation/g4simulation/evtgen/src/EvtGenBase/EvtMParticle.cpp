
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

#include "EvtGenBase/EvtMParticle.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"

EvtMParticle::EvtMParticle( int label, const EvtId& id )
{
    _id = id;
    _twospin = EvtSpinType::getSpin2( EvtPDL::getSpinType( id ) );
    _resonance.push_back( label );
}

EvtSpinAmp EvtMParticle::amplitude( const vector<EvtVector4R>& /*product*/ ) const
{
    vector<EvtSpinType::spintype> types( 2, getspintype() );
    EvtSpinAmp amp( types, EvtComplex( 0.0, 0.0 ) );

    for ( int i = -_twospin; i <= _twospin; i += 2 )
        amp( i, i ) = EvtComplex( 1.0, 0.0 );

    return amp;
}

EvtMNode* EvtMParticle::duplicate() const
{
    return new EvtMParticle( _resonance[0], _id );
}
