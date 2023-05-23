#pragma once

#include <GlobalVariables.C>

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>


#include <tpc_hits.h>

#include <stdio.h>
//#include <sstream>
#include <frog/FROG.h>

#include <string>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfun4all.so)
// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(rawtodst.so)
// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfun4allraw.so)

void Fun4All_ReadTPCHits(
  const int nEvents = 1,
  const string &fname = "/sphenix/lustre01/sphnxpro/rawdata/commissioning/tpc/pedestal/TPC_ebdc01_pedestal-00010305-0000.prdf",
  const string &outputFile = "./DSTFile.root",
  const string &outdir = "." )
  {

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  tpc_hits *TPC_HITS_set = new tpc_hits();

  se->registerSubsystem(TPC_HITS_set);

  gSystem->Load("libFROG");
  FROG *fr = new FROG();
  string inputFileName = fr->location(fname);
  cout << "Next file:" << inputFileName << endl;
  // this (DST) input manager just drives the event loop
  Fun4AllPrdfInputManager *in = new Fun4AllPrdfInputManager("PRDF1");
  in->AddFile(inputFileName);
  //in->AddFile("/sphenix/lustre01/sphnxpro/rawdata/commissioning/tpc/pedestal/TPC_ebdc01_pedestal-00010305-0000.prdf");
  se->registerInputManager(in);


  Enable::DSTOUT = true;
  Enable::DSTOUT_COMPRESS = false;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;

  if (Enable::DSTOUT)
  {
    string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
  }
  if (nEvents <= 0)
  {
    return;
  }
  cout << endl << "Running over " << nEvents << " Events" << endl;
  se->run(nEvents);
  //}
  cout << endl << "Calling End in Fun4All_ReadTPCHits.C" << endl;
  se->End();

  cout << endl << "All done, calling delete Fun4AllServer" << endl;
  delete se;

  cout << endl << "gSystem->Exit(0)" << endl;
  gSystem->Exit(0);


}