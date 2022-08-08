
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

#ifndef __EVTMPARTICLE_HH__
#define __EVTMPARTICLE_HH__

#include "EvtGenBase/EvtMNode.hh"

class EvtMParticle : public EvtMNode {
  public:
    EvtMParticle( int label, const EvtId& id );
    ~EvtMParticle() {}
    EvtSpinAmp amplitude( const vector<EvtVector4R>& product ) const override;
    int getnchild() const override { return 0; }

    EvtComplex line( const vector<EvtVector4R>& /*product*/ ) const override
    {
        return EvtComplex( 1.0, 0.0 );
    }

    EvtMNode* duplicate() const override;
};

#endif
