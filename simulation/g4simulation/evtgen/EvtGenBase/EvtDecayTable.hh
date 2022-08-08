
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

#ifndef EVTDECAYTABLE_HH
#define EVTDECAYTABLE_HH

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticleDecayList.hh"

#include <vector>

class EvtId;

// Description: Class to read in and handle the decays available
//              to EvtGen for each particle, and the model to be
//              used for each one.

class EvtDecayTable {
  public:
    static EvtDecayTable* getInstance();

    int getNMode( int ipar );

    EvtDecayBase* getDecay( int ipar, int imode );

    void readDecayFile( const std::string dec_name, bool verbose = true );
    void readXMLDecayFile( const std::string dec_name, bool verbose = true );

    bool stringToBoolean( std::string valStr );
    void checkParticle( std::string particle );

    int findChannel( EvtId parent, std::string model, int ndaug, EvtId* daugs,
                     int narg, std::string* args );

    int inChannelList( EvtId parent, int ndaug, EvtId* daugs );

    EvtDecayBase* getDecayFunc( EvtParticle* p );

    void printSummary();

    void checkConj();

    std::vector<EvtParticleDecayList> getDecayTable() { return _decaytable; };

    EvtDecayBase* findDecayModel( int aliasInt, int modeInt );
    EvtDecayBase* findDecayModel( EvtId id, int modeInt );

    bool hasPythia( int aliasInt );
    bool hasPythia( EvtId id );

    int getNModes( int aliasInt );
    int getNModes( EvtId id );

    std::vector<std::string> splitString( std::string& theString,
                                          std::string& splitter );

  protected:
    EvtDecayTable();
    ~EvtDecayTable();

  private:
    std::vector<EvtParticleDecayList> _decaytable;

    EvtDecayTable( const EvtDecayTable& ){};
    //EvtDecayTable& operator=(const EvtDecayTable&) {};
};

#endif
