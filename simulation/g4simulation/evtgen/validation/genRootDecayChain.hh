
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

#ifndef GENROOT_DECAYCHAIN_HH
#define GENROOT_DECAYCHAIN_HH

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <string>

class EvtParticle;

class genRootDecayChain {
  public:
    genRootDecayChain( const std::string& decayFileName,
                       const std::string& rootFileName,
                       const std::string& parentName, int nEvents,
                       bool storeMtmXYZ = false );

    ~genRootDecayChain();

    void run();

  protected:
    void initTree();
    void generateEvents();
    void storeDaughterInfo( EvtParticle* theParticle );
    void storeTreeInfo( EvtParticle* theParticle );
    void writeTree();

  private:
    std::string _decayFileName;
    std::string _rootFileName;
    std::string _parentName;
    int _nEvents;
    bool _storeMtmXYZ;

    TFile* _theFile;
    TTree* _theTree;
    TH1D* _probHist;
    TCanvas* _theCanvas;

    int _eventId;
    int _PDGId;
    int _dVtx;
    int _pVtx;
    int _daug;
    double _p;
    double _E;
    double _pL;
    double _EL;
    double _m;
    double _px;
    double _py;
    double _pz;
    double _pxL;
    double _pyL;
    double _pzL;
    double _t;

    int _vertexNo;
};

#endif
