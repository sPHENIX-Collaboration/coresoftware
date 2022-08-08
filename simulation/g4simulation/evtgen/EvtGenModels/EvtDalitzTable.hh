
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

#ifndef EVTDALITZTABLE_HPP
#define EVTDALITZTABLE_HPP

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include "EvtGenModels/EvtDalitzDecayInfo.hh"

#include <map>
#include <string>
#include <vector>

// Description: Model to describe a generic dalitz decay

class EvtDalitzTable {
  public:
    static EvtDalitzTable* getInstance( const std::string dec_name = "",
                                        bool verbose = true );

    bool fileHasBeenRead( const std::string dec_name );
    void readXMLDecayFile( const std::string dec_name, bool verbose = true );
    void checkParticle( std::string particle );

    void addDecay( EvtId parent, const EvtDalitzDecayInfo& dec );
    void copyDecay( EvtId parent, EvtId* daughters, EvtId copy, EvtId* copyd );

    std::vector<EvtDalitzDecayInfo> getDalitzTable( const EvtId& parent );

  protected:
    EvtDalitzTable();
    ~EvtDalitzTable();

  private:
    EvtDalitzReso getResonance( std::string shape, EvtDalitzPlot dp,
                                EvtCyclic3::Pair angPair,
                                EvtCyclic3::Pair resPair,
                                EvtSpinType::spintype spinType, double mass,
                                double width, double FFp, double FFr,
                                double alpha, double aLass, double rLass,
                                double BLass, double phiBLass, double RLass,
                                double phiRLass, double cutoffLass );
    int getDaughterPairs(
        EvtId* resDaughter, EvtId* daughter,
        std::vector<std::pair<EvtCyclic3::Pair, EvtCyclic3::Pair>>& angAndResPairs );

    std::map<EvtId, std::vector<EvtDalitzDecayInfo>> _dalitztable;
    std::vector<std::string> _readFiles;

    EvtDalitzTable( const EvtDalitzTable& );
    EvtDalitzTable& operator=( const EvtDalitzTable& );

    //to calculate probMax
    double calcProbMax( EvtDalitzPlot dp, EvtDalitzDecayInfo* model );
    double calcProb( EvtDalitzPoint point, EvtDalitzDecayInfo* model );
};

#endif
