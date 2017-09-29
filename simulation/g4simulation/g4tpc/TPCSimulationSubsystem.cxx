#include <iostream>
#include <TH1F.h>
#include <TStopwatch.h>
#include <g4detectors/PHG4Cellv1.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>
#include <TPCPadMap.h>
#include <TPCHitsContainer.h>
#include <TPCDigitsContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>

#include "TPCSimulationSubsystem.h"
#include "TPCSimulation.h"

TPCSimulationSubsystem::TPCSimulationSubsystem() :
  SubsysReco("TPCSimulationSubsystem"),
  fSimulation(new TPCSimulation()),
  fStopwatch(NULL),
  fHist(NULL),
  fTreeFileName("")
{
}

TPCSimulationSubsystem::~TPCSimulationSubsystem()
{
  if(fStopwatch) delete fStopwatch;
  delete fSimulation;
}
//======
int TPCSimulationSubsystem::Init(PHCompositeNode* node)
{
  if(verbosity>0) {
    fStopwatch = new TStopwatch();
    Fun4AllServer *se = Fun4AllServer::instance();
    fHist = new TH1F("TPCSIMTIME","TPCSIMTIME;sec/event",100,0,500);
    se->registerHisto( fHist ); // fHist now owned by me anymore
  }
  if(fTreeFileName.Length()>0) {
    fSimulation->PrepareTree(fTreeFileName);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//======
int TPCSimulationSubsystem::InitRun(PHCompositeNode *node)
{
  fSimulation->SetVerbosity(verbosity);
  PHNodeIterator iter(node);
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if(!dstNode){
    std::cout << PHWHERE << "ERROR!" << std::endl;
    exit(1);
  }
  TPCHitsContainer *hits = findNode::getClass<TPCHitsContainer>( node, "TPCHits" );
  if(!hits) {
    std::cout << PHWHERE << "ERROR!" << std::endl;
    exit(1);
  }
  fSimulation->ConnectHits( hits );

  TPCDigitsContainer *digits = findNode::getClass<TPCDigitsContainer>( node , "TPCDigits" );
  if( !digits ) {
    digits = new TPCDigitsContainer();
    dstNode->addNode( new PHIODataNode<TPCDigitsContainer>( digits, "TPCDigits", "TPCDigitsContainer" ));
  }
  fSimulation->ConnectDigits( digits );

  // creating sphenix storage element and connecting to it
  TPCPadMap *map = new TPCPadMap();
  PHG4CellContainer *phcells =  findNode::getClass<PHG4CellContainer>( node , "G4CELL_SVTX" );
  if( !phcells ) {
    phcells = new PHG4CellContainer();
    dstNode->addNode( new PHIODataNode<PHG4CellContainer>( phcells, "G4CELL_SVTX", "PHG4CellContainer" ));
  }
  for(Module_t mod=0; mod!=map->GetNumberOfModules(); ++mod) {
    for(Pad_t pad=0; pad!=map->GetNumberOfPads(mod); ++pad) {
      float rads = (map->GetRP(mod,pad)).first;
      unsigned int layers = kNPadRowsPerModule * kNSections;
      float steplyr = (kGasOuterRadius-kGasOuterRadius)/int(layers);
      unsigned int lyr = (rads-kGasInnerRadius)/steplyr;
      PHG4CellDefs::keytype key = PHG4CellDefs::SizeBinning::genkey(static_cast<unsigned short> (lyr),
								    static_cast<unsigned short> (mod),
								    static_cast<unsigned short> (pad) );
      PHG4Cellv1 *cell = new PHG4Cellv1(key);
      phcells->AddCell( cell );
      TPCDigit *dig = new TPCDigit();
      dig->SetCell( cell );
      digits->AddOrSetDigit(mod,pad,dig);
    }
    std::cout << "Pushed " << digits->GetNDigits(mod) << " in module " << int(mod) << std::endl;
  }
  delete map;
  return Fun4AllReturnCodes::EVENT_OK;
}
//======
int TPCSimulationSubsystem::process_event(PHCompositeNode *topNode)
{
  if(verbosity>1) {
    fStopwatch->Reset();
    fStopwatch->Start();
  }
  //-----
  fSimulation->Hits2Digits();
  //-----
  if(verbosity>0) {
    fHist->Fill( fStopwatch->RealTime() );
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//======
int TPCSimulationSubsystem::End(PHCompositeNode *topNode)
{
  if(fTreeFileName.Length()>0) {
    fSimulation->WriteFile();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


