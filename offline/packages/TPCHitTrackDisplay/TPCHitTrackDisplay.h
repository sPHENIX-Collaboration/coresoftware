

#ifndef __TPCHitTrackDisplay_H__
#define __TPCHitTrackDisplay_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>
#include <fstream>


//Forward declerations
class PHCompositeNode;
class TFile; 
class TTree;
class RawTowerContainer;
class RawTowerGeomContainer;
class TH2F;
class TProfile;
class SvtxTrackMap;


RawTowerContainer *cemctowers;
RawTowerContainer *hcalotowers;
RawTowerContainer *hcalitowers;

RawTowerGeomContainer *cemctowergeom;
RawTowerGeomContainer *hcalotowergeom;
RawTowerGeomContainer *hcalitowergeom;

// Writes json file to be used to display an event with:
// https://www.sphenix.bnl.gov/edisplay/
class TPCHitTrackDisplay: public SubsysReco
{
 public: 
  //Default constructor
  TPCHitTrackDisplay(const std::string &name="TPCHitTrackDisplay"/*, bool &tpcRaw=True*/);

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

  int InitRun(PHCompositeNode *); 

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int EndRun(const int);
  
  void set_pdgcode(const int thispdgcode) { _pdgcode = thispdgcode; }

  // set the ADC cut for displaying trackless clusters
  void setCutADC(float value) { m_cut_ADC = value; }

  // Boolean for whether or not to include clusters without an associted track above a certain ADC value 
  void setIncludeTracklessClusters(float value) { m_trackless_clusters = value; }
   
 private:

  float m_cut_ADC;  
  bool m_trackless_clusters; 

  //Event counter
  int _event;
  int _pdgcode;
  std::string _fileName;  
 
  //bool isRawData;
    
  //User modules
  void SimulationOut(PHCompositeNode*);
  //void TPCRawOut(PHCompositeNode*);

};

#endif //* __TPCHitTrackDisplay_H__ *//
