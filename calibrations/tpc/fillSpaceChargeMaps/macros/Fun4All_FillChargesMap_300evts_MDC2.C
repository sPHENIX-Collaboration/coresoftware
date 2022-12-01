#pragma once
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <frog/FROG.h>

#include <fillSpaceChargeMaps.h>

#include <stdio.h>
//#include <sstream>

#include <string>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfillSpaceChargeMaps.so)
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
void Fun4All_FillChargesMap_300evts_MDC2(  const int nEvents = 10, const int eventsInFileStart = 0, const int eventsBeamCrossing = 1508071, const string &fname = "/sphenix/sim/sim01/sphnxpro/Micromegas/2/G4Hits_sHijing_0-12fm_000000_001000.root", const string &foutputname = "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_000000_001000.root" )
{
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////


  //char fname[100];
  //char foutputname[100];
  //for (int i=0;i<nFiles;i++){
  //int eventsInFileStart = i*1000;
  //int eventsInFileEnd = (i+1)*1000;
  //sprintf(fname, "/sphenix/sim/sim01/sphnxpro/Micromegas/2/G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd);
  //sprintf(foutputname, "/sphenix/user/shulga/Work/IBF/DistortionMap/Files/slim_G4Hits_sHijing_0-12fm_%06d_%06d.root",eventsInFileStart,eventsInFileEnd);
  std::vector<int> bXs = readBeamXings();
  std::vector<int> bXs_sel;
  //std::vector<int> bXs_sel_end;
  std::vector<int>::iterator it = std::find(bXs.begin(), bXs.end(), eventsBeamCrossing);
  int index=0; 
  index = std::distance(bXs.begin(), it);
  cout<<"Index="<<index<<endl;
  for(int n=0;n<30;n++){
    int bXN=index+n*300;
    bXs_sel.push_back(bXs[bXN]);
    //int n_bX = closest(bXs, bXs[bXN]);
    //bXs_sel_end.push_back(n_bX); 
    //int id_bX = std::distance(bXs.begin(), it);
    //std::vector<int>::iterator it_bX = std::find(bXs.begin(), bXs.end(), n_bX);
    cout<<"bX="<<bXs[bXN]<<endl;
    //cout<<"bX="<<bXs[bXN]<<"last event:"<<n_bX<<" evt:"<<std::distance(bXs.begin(), it_bX)<<endl;
  }
  Fun4AllServer *se = Fun4AllServer::instance();
  string cd_name = "fillSpaceChargeMaps"+std::to_string(eventsInFileStart);
  //cout<<fname_tmp<<endl;
  fillSpaceChargeMaps *dist_calc = new fillSpaceChargeMaps(cd_name, foutputname);
  dist_calc->SetFrequency(50);
  dist_calc->SetEvtStart(eventsInFileStart);
  //dist_calc->SetBeamXing(eventsBeamCrossing); // Set beam crosssing bias
  dist_calc->SetBeamXing(bXs_sel);// Set beam crosssing bias
  //dist_calc->SetBeamXingEnd(bXs_sel_end);// Set last beam crosssing for the biases
  //dist_calc->SetAvg(1); //Set average calculation
  dist_calc->SetUseIBFMap(false);//false);
  //dist_calc->SetGain(2e3*48.7/71.5);
  dist_calc->SetGain(1400);
  dist_calc->SetIBF(0.004);
  dist_calc->UseSliming(0);//Turn off TTree filling and recording

  dist_calc->UseFieldMaps(1);//1); //setting field maps to shift electron position
  //Set pp colliding system
  dist_calc->SetCollSyst(0); //setting pp with = 1

  se->registerSubsystem(dist_calc);
  

  gSystem->Load("libFROG");
  FROG *fr = new FROG();
  string inputFileName = fr->location(fname);
  cout << "Next file:" << inputFileName << endl;
  // this (DST) input manager just drives the event loop
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(inputFileName);
  se->registerInputManager(in);
  // events = 0 => run till end of input file
  if (nEvents <= 0)
  {
    return;
  }
  cout << endl << "Running over " << nEvents << " Events" << endl;
  se->run(nEvents);

  cout << endl << "Calling End in Fun4All_fillSpaceChargeMaps.C" << endl;
  se->End();

  cout << endl << "All done, calling delete Fun4AllServer" << endl;
  delete se;

  cout << endl << "gSystem->Exit(0)" << endl;
  gSystem->Exit(0);
}