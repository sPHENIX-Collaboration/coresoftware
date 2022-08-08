
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

#ifndef __EVTMTREE_HH__
#define __EVTMTREE_HH__

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtMNode.hh"
#include "EvtGenBase/EvtMParticle.hh"
#include "EvtGenBase/EvtMRes.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSpinAmp.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <vector>
using std::vector;

#include <string>
using std::string;

typedef string::const_iterator ptype;

class EvtParticle;

class EvtMTree {
  public:
    EvtMTree( const EvtId*, unsigned int );
    ~EvtMTree();

    // return the invariant amplitude of the entire tree
    EvtSpinAmp amplitude( EvtParticle* ) const;

    // add a decay tree to the list of trees that we posess
    void addtree( const string& );

  private:
    vector<EvtMNode*> _root;
    vector<string> _lbltbl;
    double _norm;

    bool parsecheck( char, const string& );
    void parseerror( bool, ptype&, ptype&, ptype& );

    string parseId( ptype&, ptype&, ptype& );
    string parseKey( ptype&, ptype&, ptype& );
    vector<string> parseArg( ptype&, ptype&, ptype& );
    vector<EvtComplex> parseAmps( ptype&, ptype&, ptype& );
    vector<EvtMNode*> duplicate( const vector<EvtMNode*>& ) const;
    vector<vector<EvtMNode*>> unionChildren( const string&,
                                             vector<vector<EvtMNode*>>& );
    vector<vector<EvtMNode*>> parseChildren( ptype&, ptype&, ptype& );
    vector<EvtMNode*> parsenode( const string&, bool );
    bool validTree( const EvtMNode* ) const;

    vector<EvtMNode*> makeparticles( const string& );
    EvtMRes* makeresonance( const EvtId&, const string&, const vector<string>&,
                            const string&, const vector<EvtComplex>&,
                            const vector<EvtMNode*>& );

    EvtSpinAmp getrotation( EvtParticle* ) const;
};

#endif
