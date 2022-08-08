
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

#ifndef __EVTMRES_HH__
#define __EVTMRES_HH__

#include "EvtGenBase/EvtMNode.hh"

class EvtMRes;

class EvtMLineShape {
  public:
    virtual ~EvtMLineShape(){};

    void setres( EvtMRes* n ) { _node = n; }
    virtual EvtComplex shape( const vector<EvtVector4R>& product ) const = 0;

    virtual EvtMLineShape* duplicate() const = 0;

  protected:
    EvtMRes* _node;
};

class EvtMRes : public EvtMNode {
  public:
    ~EvtMRes();

    int getnchild() const override { return _children.size(); }

    EvtComplex line( const vector<EvtVector4R>& product ) const override
    {
        return _lineshape->shape( product );
    }

  protected:
    // store the child nodes
    vector<EvtMNode*> _children;

    // store the parametrization amplitudes in some kind
    EvtSpinAmp _amp;

    // store the lineshape of the resonance
    EvtMLineShape* _lineshape;
};

#endif
