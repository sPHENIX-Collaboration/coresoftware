#include "CaloVtxReco.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jet.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfoDefs.h>
#include <globalvertex/CaloVertexv1.h>
#include <globalvertex/CaloVertexMapv1.h>

#include <cmath>
static const float radius_EM = 93.5;
static const float radius_OH = 225.87;

//____________________________________________________________________________..
CaloVtxReco::CaloVtxReco(const std::string &name, const std::string &jetnodename, const int debug):
  SubsysReco(name), _name(name), _jetnodename(jetnodename), _debug(debug)
{  
  _calovtxmap = NULL;
}

//____________________________________________________________________________..
CaloVtxReco::~CaloVtxReco()
= default;

int CaloVtxReco::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "RUN Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator dstiter(dstNode);
  
  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode)
  {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  _calovtxmap = findNode::getClass<CaloVertexMap>(globalNode, "CaloVertexMap");
  if (!_calovtxmap)
  {
    _calovtxmap = new CaloVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(_calovtxmap, "CaloVertexMap", "PHObject");
    globalNode->addNode(VertexMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloVtxReco::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) { std::cout << "Initializing!" << std::endl;
}
  if(createNodes(topNode) == Fun4AllReturnCodes::ABORTRUN)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}


float CaloVtxReco::new_eta(int channel, TowerInfoContainer* towers, RawTowerGeomContainer* geom, RawTowerDefs::CalorimeterId caloID, float testz)
{
  int key = towers->encode_key(channel);
  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(caloID, towers->getTowerEtaBin(key), towers->getTowerPhiBin(key));

  RawTowerGeom* tower_geom = geom->get_tower_geometry(geomkey);
  float oldeta = tower_geom->get_eta();
  
  float radius = (caloID==RawTowerDefs::CalorimeterId::HCALIN?radius_EM:radius_OH);
  float towerz = radius/(tan(2*atan(std::exp(oldeta))));
  float newz = towerz + testz;
  float newTheta = std::atan2(radius,newz);
  float neweta = -log(tan(0.5*newTheta));
  
  return neweta;
}

float get_dphi(float phi1, float phi2)
{
  float dphi = abs(phi1-phi2);
  if(dphi > M_PI) { dphi = 2*M_PI - dphi;
}
  return dphi;
}

int CaloVtxReco::process_event(PHCompositeNode *topNode)
{

  

  if(_debug > 1) { std::cout << std::endl << std::endl << std::endl << "CaloVtxReco: Beginning event processing" << std::endl;
}
  _zvtx = std::numeric_limits<float>::quiet_NaN();
  JetContainer *jetcon = findNode::getClass<JetContainerv1>(topNode, _jetnodename);
  TowerInfoContainer *towers[3];
  towers[0] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  towers[1] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  towers[2] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  RawTowerGeomContainer *geom[3];
  geom[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  const int nz = 601;
  const int njet = 2;
  Jet* jets[njet];
  float jpt[njet] = {0};
  float jemsum[njet] = {0};
  float johsum[njet] = {0};
  float jemeta[njet] = {0};
  float joheta[njet] = {0};

  if(jetcon)
    {
      int tocheck = jetcon->size();
      if(_debug > 2) { std::cout << "Found " << tocheck << " jets to check..." << std::endl;
}
      for(int i=0; i<tocheck; ++i)
        {
          Jet *jet = jetcon->get_jet(i);
          if(jet)
            {
	      float pt = jet->get_pt();
	      if(pt < _jet_threshold) { continue;
}
	      if(pt > jpt[0])
		{
		  jpt[1] = jpt[0];
		  jets[1] = jets[0];
		  jets[0] = jet;
		  jpt[0] = pt;
		}
	      else if(pt > jpt[1])
		{
		  jets[1] = jet;
		  jpt[1] = pt;
		}
	    }
	}      
    }
  else
    {
      if(_debug > 0) { std::cout << "no jets" << std::endl;
}
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  if(jpt[0] == 0)
    {
      if(_debug > 2) { std::cout << "NO JETS > 5 GeV!" << std::endl;
}
    }
	

  float metric = FLT_MAX;
  for(int i=0; i<nz; ++i)
    {
      float testz = -300+i;
      float testmetric = 0;
      for(int j=0; j<njet; ++j)
	{
	  if(jpt[j] == 0) { continue;
}
	  jemsum[j] = 0;
	  johsum[j] = 0;
	  jemeta[j] = 0;
	  joheta[j] = 0;
	  for(auto comp: jets[j]->get_comp_vec())
	    {
	      if(comp.first==5 || comp.first == 26) continue;
	      unsigned int channel = comp.second;
	      if(comp.first==7 || comp.first == 27)
		{
		  TowerInfo* tower = towers[2]->get_tower_at_channel(channel);
		  if(tower->get_energy() < 0.1) continue;
		  johsum[j] += tower->get_energy();
		  float neweta = new_eta(channel, towers[2], geom[2], RawTowerDefs::CalorimeterId::HCALOUT, testz);
		  joheta[j] += neweta*tower->get_energy();
		}
	      if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		{
		  TowerInfo* tower = towers[0]->get_tower_at_channel(channel);
		  if(tower->get_energy() < 0.1) continue;
		  jemsum[j] += tower->get_energy();
		  float neweta = new_eta(channel, towers[0], geom[1], RawTowerDefs::CalorimeterId::HCALIN, testz);
		  jemeta[j] += neweta*tower->get_energy();
		}
	    }
	  jemeta[j] /= jemsum[j];
	  joheta[j] /= johsum[j];
	  if((jemsum[j] == 0 || johsum[j] == 0) && _debug > 1) { std::cout << "zero E sum in at least one calo for a jet" << std::endl;
}
	  if(!std::isnan(jemeta[j]) && !std::isnan(joheta[j]))
	    {
	      testmetric += pow(jemeta[j]-joheta[j],2);
	    }
	}
      if(_debug > 3) { std::cout << "metric: " << testmetric << std::endl;
}
      if(testmetric < metric && testmetric != 0)
	{
	  metric = testmetric;
	  _zvtx = testz;
	}
    }
  if(abs(_zvtx) < 305)
    {
      if(_debug > 2) { std::cout << "optimal z: " << _zvtx << std::endl;
}
      CaloVertex *vertex = new CaloVertexv1();
      _zvtx *= _calib_factor; //calibration factor from simulation
      vertex->set_z(_zvtx);
      _calovtxmap->insert(vertex);
      if(_debug > 3) { std::cout << "CaloVtxReco: end event" << std::endl;
}
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloVtxReco::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "CaloVtxReco::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void CaloVtxReco::Print(const std::string &what) const
{
  std::cout << "CaloVtxReco::Print(const std::string &what) const Printing info for " << what << std::endl;
}
