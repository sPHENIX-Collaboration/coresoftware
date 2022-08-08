
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

#ifndef __EVTMNODE_HH__
#define __EVTMNODE_HH__

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSpinAmp.hh"
#include "EvtGenBase/EvtSymTable.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <vector>
using std::vector;

#include <string>
using std::string;

class EvtMNode {
  public:
    EvtMNode() {}
    virtual ~EvtMNode(){};

    // calculate the amplitude associated event this->children return a
    // vector of the form A_{\lambda this} and sum over allowed angular
    // momenta of the children
    virtual EvtSpinAmp amplitude( const vector<EvtVector4R>& product ) const = 0;

    // get the 4 vector associated with this node
    EvtVector4R get4vector( const vector<EvtVector4R>& product ) const;

    // get twice the spin of the particle
    int getspin() const { return _twospin; }
    EvtSpinType::spintype getspintype() const
    {
        return EvtPDL::getSpinType( _id );
    }

    // get the id of this node
    EvtId getid() const { return _id; }

    // return which particles this is a combination of
    const vector<int>& getresonance() const { return _resonance; }

    void setparent( EvtMNode* parent ) { _parent = parent; }
    EvtMNode* getparent() const { return _parent; }

    // get the number of children that this node has
    virtual int getnchild() const = 0;

    // return the value of the resonance shape
    virtual EvtComplex line( const vector<EvtVector4R>& product ) const = 0;

    // return a pointer node
    virtual EvtMNode* duplicate() const = 0;

  protected:
    // store the EvtId of the particle (just in case we need it to access
    // further informatoin about it)
    EvtId _id;

    // store TWICE the spin of this resonance (this is to deal with spin 1/2
    int _twospin;

    // store the particles that form this resonance, this should match up
    // with the child nodes from below, and is calculated internally
    vector<int> _resonance;

    // store the parent node of this one
    EvtMNode* _parent;
};

#endif
