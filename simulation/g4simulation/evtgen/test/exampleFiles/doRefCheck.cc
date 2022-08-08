
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

//
// Script to automate reference comparisons
//

#include "TCanvas.h"
#include "TCollection.h"
#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"
#include "TROOT.h"
#include "TStyle.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

void doRefCheck( TString file, bool makepsfiles = false, bool writeRef = false )
{
    // read in root file
    TFile* f1 = new TFile( file, "READ" );
    //  f1->ls();
    TKey* key;
    TIter nextkey( f1->GetListOfKeys() );

    std::vector<std::string> names;

    while ( key = (TKey*)nextkey() ) {
        if ( strcmp( key->GetClassName(), "TH1F" ) != 0 )
            continue;
        TH1F* h1 = (TH1F*)key->ReadObj();
        names.push_back( ( h1->GetName() ) );
    }

    if ( names.size() == 0 ) {
        std::cout << "NO histograms found?? Missing file?\n";
        return;
    }

    // read in references if they are found.
    std::map<std::string, TH1F*> references;
    TString refFile = file;
    refFile += ".ref";
    ifstream* infile = new ifstream( refFile.Data(), std::ios::in );
    char tChar[100];
    std::vector<double> datums;
    std::string refTitle;
    double refLow = 0.;
    double refHigh = 0.;

    while ( infile->good() ) {
        infile->getline( tChar, 100 );
        if ( infile->good() ) {
            if ( strcmp( tChar, "start histogram" ) == 0 ) {
                if ( datums.size() > 0 ) {
                    TString tTitle = "ref";
                    tTitle += references.size();
                    TH1F* hTemp = new TH1F( tTitle, tTitle, datums.size() - 2,
                                            refLow, refHigh );
                    for ( unsigned int i = 0; i < datums.size(); i++ ) {
                        hTemp->SetBinContent( i, datums[i] );
                        //std::cout << datums[i] << std::endl;
                    }
                    references[refTitle] = hTemp;
                    datums.clear();
                }
                infile->getline( tChar, 100 );
                refTitle = std::string( tChar );
                infile->getline( tChar, 100 );
                refLow = atof( tChar );
                infile->getline( tChar, 100 );
                refHigh = atof( tChar );
            } else {
                datums.push_back( atof( tChar ) );
            }
        }
    }
    if ( datums.size() > 0 ) {
        TString tTitle = "ref";
        tTitle += references.size();
        TH1F* hTemp = new TH1F( tTitle, tTitle, datums.size() - 2, refLow,
                                refHigh );
        for ( unsigned int i = 0; i < datums.size(); i++ )
            hTemp->SetBinContent( i, datums[i] );
        references[refTitle] = hTemp;
    }

    gROOT->SetStyle( "Plain" );    // make plot background truly white
    gStyle->SetPalette( 1 );    // palette set from violet to red, root chooses scale

    fstream* outf = 0;
    if ( writeRef ) {
        TString refOut = file;
        refOut += ".ref.new";
        outf = new fstream( refOut.Data(), std::ios::out );
    }

    for ( unsigned int i = 0; i < names.size(); i++ ) {
        TString canName = "can";
        canName += i;
        TCanvas* can = new TCanvas( canName, canName, 600, 400 );

        TH1F* h1 = (TH1F*)f1->Get( names[i].c_str() );
        h1->SetMinimum( 0. );
        h1->Draw();

        TH1F* h2 = references[std::string( h1->GetTitle() )];
        std::cout << "Histogram: " << h1->GetTitle();
        if ( h2 != 0 ) {
            h2->SetLineColor( kRed );
            h2->Draw( "same" );

            if ( h1->GetNbinsX() != h2->GetNbinsX() ) {
                std::cout << " reference has different binning! "
                          << h1->GetNbinsX() << " " << h2->GetNbinsX()
                          << std::endl;
            } else {
                int nDiffs = 0;
                for ( int i = 0; i < h1->GetNbinsX(); i++ )
                    if ( h1->GetBinContent( i ) != h2->GetBinContent( i ) )
                        nDiffs++;
                if ( nDiffs > 0 )
                    std::cout << " has " << nDiffs
                              << " bins that were different from reference"
                              << std::endl;
                else
                    std::cout << " agrees with reference\n";
            }
        } else {
            std::cout << " no reference found\n";
        }

        TString title = h1->GetTitle();
        if ( outf != 0 ) {
            ( *outf ) << "start histogram\n";
            ( *outf ) << title.Data() << std::endl;
            ( *outf ) << h1->GetBinLowEdge( 1 ) << std::endl
                      << h1->GetBinLowEdge( 1 + h1->GetNbinsX() ) << std::endl;
            for ( int j = 0; j <= h1->GetNbinsX() + 1; j++ )
                ( *outf ) << h1->GetBinContent( j ) << std::endl;
        }

        if ( makepsfiles ) {
            TString psName = file;
            psName += "_";
            psName += i;
            psName += ".ps";
            can->Print( psName );
        }
    }

    if ( outf != 0 )
        outf->close();
}
