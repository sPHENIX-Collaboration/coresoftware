
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

#include "EvtGenBase/EvtMBreitWigner.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>

using std::endl;

EvtMBreitWigner::EvtMBreitWigner( const EvtId& id, const vector<string>& args )
{
    if ( args.size() != 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Unknown input arguments passed in to lineshape." << endl;
        ::abort();
    }

    _id = id;
    _width = EvtPDL::getWidth( id );
    _resmass = EvtPDL::getMeanMass( id );
}

EvtComplex EvtMBreitWigner::shape( const vector<EvtVector4R>& product ) const
{
    static EvtComplex I( 0.0, 1.0 );
    double mass = _node->get4vector( product ).mass();

    return sqrt( _width / ( EvtConst::twoPi ) ) * 1 /
           ( mass - _resmass - I * _width / 2 );
}

EvtMLineShape* EvtMBreitWigner::duplicate() const
{
    vector<string> args;
    EvtMLineShape* tmp = new EvtMBreitWigner( _id, args );
    return tmp;
}
