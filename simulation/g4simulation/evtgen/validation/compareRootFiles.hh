
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

#ifndef COMPARE_ROOTFILES_HH
#define COMPARE_ROOTFILES_HH

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

class compareRootFiles {
  public:
    compareRootFiles( string decayName, string newFileName, string oldFileName,
                      string description );
    ~compareRootFiles();

    void getNDaugPlots();
    void getAllIdPlots();
    void getPartIdPlots();
    void getMtmPlots();

  private:
    TH1D* getDaugHist( TFile* theFile, string histName );
    TH1D* getAllIdHist( TFile* theFile, string histName );
    TH1D* getPartIdHist( TFile* theFile, string histName );
    TH1D* getMtmHist( TFile* theFile, string histName, vector<int> groupInts );

    TH1D* _emptyHist;

    int getPartGroup( int PDGId );

    string _decayName, _newFileName, _oldFileName, _description;

    TFile* _newFile;
    TFile* _oldFile;

    TCanvas* _theCanvas;
    int _nGroupMax;
};

#endif
