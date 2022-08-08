
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

#ifndef EVTPDL_HH
#define EVTPDL_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPartProp.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtStringHash.hh"

#include <map>
#include <vector>
#include <string>
#include <iostream>

const int SPIN_NAME_LENGTH = 100;

class EvtPDL final {
  public:
    EvtPDL();

    void read( const std::string& fname );
    void readPDT( std::istream& data );

    static double getMeanMass( EvtId i );
    static double getMass( EvtId i );
    static double getRandMass( EvtId i, EvtId* parId, int nDaug, EvtId* dauId,
                               EvtId* othDaugId, double maxMass,
                               double* dauMasses );
    static double getMassProb( EvtId i, double mass, double massPar, int nDaug,
                               double* massDau );

    static double getMaxMass( EvtId i );
    static double getMinMass( EvtId i );
    //the number we got from PDT
    static double getMaxRange( EvtId i );
    static double getWidth( EvtId i );
    static double getctau( EvtId i );
    static int getStdHep( EvtId id );
    static int getLundKC( EvtId id );

    // Function to retrieve EvtId from PythiaID
    static EvtId evtIdFromLundKC( int pythiaId );
    static EvtId evtIdFromStdHep( int stdhep );
    static EvtId chargeConj( EvtId id );
    static int chg3( EvtId i );
    static EvtSpinType::spintype getSpinType( EvtId i );
    static EvtId getId( const std::string& name );
    static std::string name( EvtId i );
    static void alias( EvtId num, const std::string& newname );
    static void aliasChgConj( EvtId a, EvtId abar );
    static size_t entries();
    static EvtId getEntry( int i );
    static void reSetMass( EvtId i, double mass );
    static void reSetWidth( EvtId i, double width );
    static void reSetMassMin( EvtId i, double mass );
    static void reSetMassMax( EvtId i, double mass );
    static void reSetBlatt( EvtId i, double blatt );
    static void reSetBlattBirth( EvtId i, double blatt );
    static void includeBirthFactor( EvtId i, bool yesno );
    static void includeDecayFactor( EvtId i, bool yesno );
    static void changeLS( EvtId i, std::string& newLS );
    static void setPWForDecay( EvtId i, int spin, EvtId d1, EvtId d2 );
    static void setPWForBirthL( EvtId i, int spin, EvtId par, EvtId othD );

  private:
    void setUpConstsPdt();

    static unsigned int _firstAlias;
    static int _nentries;

    static std::vector<EvtPartProp>& partlist()
    {
        static std::vector<EvtPartProp> s_partlist;
        return s_partlist;
    }

    static std::map<std::string, int> _particleNameLookup;

};    // EvtPDL.h

#endif
