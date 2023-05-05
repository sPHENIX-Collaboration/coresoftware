#pragma once
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <readDigitalCurrents.h>

#include <stdio.h>
//#include <sstream>
#include <frog/FROG.h>

#include <string>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfun4all.so)
// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libreadDigitalCurrents.so)
// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libg4dst.so)

void Fun4All_FillDCMap(  const int nEvents = 1000, const int eventsInFileStart = 0, const int eventsBeamCrossing = 1508071, const string &fname = "/sphenix/user/shulga/Work/IBF/macros/detectors/sPHENIX/Files/DST_G4Hits_sHijing_0-12fm_005000_006000.root", const string &foutputname = "./Files/hists_G4Hits_sHijing_0-12fm_000000_001000.root" )//DST_G4sPHENIX_1000evt.root")//G4sPHENIX.root" )
{
  // /sphenix/user/frawley/new_macros_april27/macros/detectors/sPHENIX/Reconstructed_DST_Hijing_50kHz_00000.root
  
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////

  Fun4AllServer *se = Fun4AllServer::instance();
  string cd_name = "readDigitalCurrents"+std::to_string(eventsInFileStart);
  //cout<<fname_tmp<<endl;
  readDigitalCurrents *dist_calc = new readDigitalCurrents(cd_name, foutputname);
  //readDigitalCurrents *dist_calc = new readDigitalCurrents();
  se->registerSubsystem(dist_calc);
  dist_calc->SetEvtStart(eventsInFileStart);
  dist_calc->SetBeamXing(eventsBeamCrossing); // Set beam crosssing bias
  dist_calc->SetCollSyst(0); //setting pp with = 1
  dist_calc->SetIBF(0.004);
  dist_calc->SetCCGC(1);//to use PHG4CylinderCellGeom

  gSystem->Load("libFROG");
  FROG *fr = new FROG();
  string inputFileName = fr->location(fname);
  cout << "Next file:" << inputFileName << endl;
  // this (DST) input manager just drives the event loop
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(inputFileName);//fname);
  se->registerInputManager(in);
  // events = 0 => run till end of input file
  if (nEvents <= 0)
  {
    return;
  }
  cout << endl << "Running over " << nEvents << " Events" << endl;
  se->run(nEvents);
  //}
  cout << endl << "Calling End in Fun4All_readDigitalCurrents.C" << endl;
  se->End();

  cout << endl << "All done, calling delete Fun4AllServer" << endl;
  delete se;

  cout << endl << "gSystem->Exit(0)" << endl;
  gSystem->Exit(0);
}
