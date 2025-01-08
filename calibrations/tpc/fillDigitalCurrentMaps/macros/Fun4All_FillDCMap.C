#pragma once
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <readDigitalCurrents.h>

#include <stdio.h>
#include <frog/FROG.h>

#include <string>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libreadDigitalCurrents.so)
R__LOAD_LIBRARY(libg4dst.so)

std::vector<int> readBeamXings();

std::vector<int> readBeamXings(){
  //cout << "fillSpaceChargeMaps::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << endl;
  std::vector<int> bXs;
  string line;
  string txt_file = "./data/timestamps_50kHz_1M.txt";
  ifstream InputFile (txt_file);
  //std::map<int,int> timestamps;
  if (InputFile.is_open()){
    int n_line=0;
    while ( getline (InputFile,line) )
    {
        n_line++;
      //cout << line << '\n';
      if(n_line>3){
        std::istringstream is( line );
        double n[2] = {0,0};
        int i = 0;
        while( is >> n[i] ) {    
            i++;    
        }
        //_timestamps[n[0]]=n[1];
        bXs.push_back(int(n[0]));
      }
    }
    InputFile.close();
  }
 return bXs;
}

int closest(std::vector<int> const& vec, int value) {
    auto const it = std::lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end()) { return -1; }

    return *it;
}

void Fun4All_FillDCMap(  const int nEvents = 1000, const int eventsInFileStart = 0, const int eventsBeamCrossing = 1508071, const string &fname = "/sphenix/user/shulga/Work/IBF/macros/detectors/sPHENIX/Files/DST_G4Hits_sHijing_0-12fm_005000_006000.root", const string &foutputname = "./Files/hists_G4Hits_sHijing_0-12fm_000000_001000.root" )//DST_G4sPHENIX_1000evt.root")//G4sPHENIX.root" )
{
  // /sphenix/user/frawley/new_macros_april27/macros/detectors/sPHENIX/Reconstructed_DST_Hijing_50kHz_00000.root
  
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  std::vector<int> bXs = readBeamXings();
  std::vector<int> bXs_sel;

  std::vector<int>::iterator it = std::find(bXs.begin(), bXs.end(), eventsBeamCrossing);
  int index=0; 
  index = std::distance(bXs.begin(), it);
  cout<<"Index="<<index<<endl;
  for(int n=0;n<30;n++){
    int bXN=index+n*300;
    bXs_sel.push_back(bXs[bXN]);
    cout<<"bX="<<bXs[bXN]<<endl;
  }

  Fun4AllServer *se = Fun4AllServer::instance();
  string cd_name = "readDigitalCurrents"+std::to_string(eventsInFileStart);
  //cout<<fname_tmp<<endl;
  readDigitalCurrents *dist_calc = new readDigitalCurrents(cd_name, foutputname);
  //readDigitalCurrents *dist_calc = new readDigitalCurrents();
  se->registerSubsystem(dist_calc);
  dist_calc->SetEvtStart(eventsInFileStart);
  dist_calc->SetBeamXing(bXs_sel);// Set beam crosssing bias
  //dist_calc->SetBeamXing(eventsBeamCrossing); // Set beam crosssing bias
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
