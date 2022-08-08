
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

#include "EvtGenBase/EvtMNode.hh"

#include "EvtGenBase/EvtPatches.hh"

EvtVector4R EvtMNode::get4vector( const vector<EvtVector4R>& product ) const
{
    EvtVector4R res( 0.0, 0.0, 0.0, 0.0 );
    vector<int>::const_iterator iter;

    for ( iter = _resonance.begin(); iter != _resonance.end(); ++iter )
        res += product[*iter];

    return res;
}
