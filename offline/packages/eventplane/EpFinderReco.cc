#include "EpFinderReco.h"
#include "EpFinder.h"
#include "EpInfo.h"    // for EpInfo
#include "EpInfov1.h"  // for EpInfo

#include <epd/EpdGeomV1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoDefs.h>

#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>

#include <cdbobjects/CDBHistos.h>

//#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH3.h>  
#include <TSystem.h>
#include <TRandom.h>
#include <TVector3.h>  // for TVector3

#include <algorithm>  // for max
#include <cmath>      // for M_PI
#include <cstdlib>    // for NULL, exit, getenv
#include <iostream>
#include <map>  // for _Rb_tree_const_iterator
#include <utility>
#include <vector>  // for vector
#include <fstream>
 
EpFinderReco::EpFinderReco(const std::string &name)
  : SubsysReco(name)
  , detector("NONE")
{
}

EpFinderReco::~EpFinderReco()
{
  for(int i = 0; i < 2; i++)
  {
    delete EpFinder_det[i];
  }
}

int EpFinderReco::Init(PHCompositeNode *topNode)
{
    
  for(int i = 0; i < 2; i++)
  {
     EpFinder_det[i] = new EpFinder(1,3,i);
     EpFinder_det[i]->Report();
  }

    //Load sEPD truncation hist if it exists
    CDBHistos *cdbhistos = new CDBHistos(truncationFile);
    cdbhistos->LoadCalibrations();
  
    if(!(cdbhistos->getHisto(truncationhist)))
    {
        mTruncationInput = NULL;
    }
    else
    {
        mTruncationInput = dynamic_cast<TH3 *> (cdbhistos->getHisto(truncationhist));
    }

   return CreateNodes(topNode);
}

int EpFinderReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *AlgoNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", _algonode));
  if (!AlgoNode)
  {
    AlgoNode = new PHCompositeNode(_algonode);
    dstNode->addNode(AlgoNode);
  }

  EpNode += detector;
  if((detector == "EPD") || (detector == "BBC"))
  {
    EpNode += "_SOUTH";
    EventPlaneNodeName.push_back(EpNode);
    EpNode = "EPINFO_" + detector;
    EpNode += "_NORTH";
    EventPlaneNodeName.push_back(EpNode);
  }
  else
  {
    EpNode = "EPINFO_" + detector;
    EventPlaneNodeName.push_back(EpNode);
  }
    
 for (unsigned int i = 0; i < EventPlaneNodeName.size(); i++)
 {
    EpInfo *EpInfo_det = new EpInfov1();
    PHIODataNode<PHObject> *EpInfo_det_node = new PHIODataNode<PHObject>(EpInfo_det, EventPlaneNodeName[i], "PHObject");
    AlgoNode->addNode(EpInfo_det_node);
 }

 return Fun4AllReturnCodes::EVENT_OK;

}

int EpFinderReco::process_event(PHCompositeNode *topNode)
{

  GetNodes(topNode);
  GetEventPlanes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int EpFinderReco::End(PHCompositeNode * /*topNode*/)
{
    for(int i = 0; i < 2; i++)
    {
      EpFinder_det[i]->Finish(i);
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int EpFinderReco::ResetEvent(PHCompositeNode * /*topNode*/)
{
    
  for(int i = 0; i < 2; i++)
  {
    EpFinder_det[i]->ResetEvent();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//

void EpFinderReco::GetEventPlanes(PHCompositeNode *topNode)
{
   
//next PR will add MBD reco, ignoring use of MBD g4hits for this iteration 
 
 if(detector != "EPD")
 {
     std::cout << detector <<"\t is not currently implemented, choose EPD for detector choice \t"<<std::endl;
     exit(1);
 }
 else
 {
    //for sEPD phi weighting
    static bool first = true;

      if(first){
           
          for(int i=0; i<16; i++){
              Nepd_phi_list[i].clear();
              Sepd_phi_list[i].clear();
          }
          
          for(int i=0; i<1; i++){
              Nepd_phi_list0[i].clear();
              Sepd_phi_list0[i].clear();
          }
          
          for(int i=1; i<16; i++){
             for(int j=0; j<24; j++){
               std::pair<int,int> newPair(i,j);
                 Nepd_phi_list[i].push_back(newPair);
                 Sepd_phi_list[i].push_back(newPair);
             }
           }
          
          for(int i=0; i<1; i++){
            for(int j=0; j<12; j++){
              std::pair<int,int> newPair(i,j);
                Nepd_phi_list0[i].push_back(newPair);
                Sepd_phi_list0[i].push_back(newPair);
            }
          }
          
          first = false;
       }
 
   CentralityInfov1 *cent = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");
   if (!cent)
   {
      std::cout << " ERROR -- can't find CentralityInfo node" << std::endl;
       exit(1);
   }

   if(m_CentType == kCentrality::kBimp) 
   {
      cent_index = cent->get_centile(CentralityInfo::PROP::bimp)/10;
   } 
   else if(m_CentType == kCentrality::kEPD)
   {
      cent_index = cent->get_centile(CentralityInfo::PROP::epd_NS)/10;
   }
   else if(m_CentType == kCentrality::kMBD)
   {
      cent_index = cent->get_centile(CentralityInfo::PROP::mbd_NS)/10;
   }

   TowerInfoContainerv1 *_towerinfos = findNode::getClass<TowerInfoContainerv1>(topNode, TowerNode + detector);
   if (!_towerinfos)
   {
      std::cout << "Could not locate tower info node " << TowerNode + detector << std::endl;
      exit(1);
   }

   EpdGeom *_epdgeom = findNode::getClass<EpdGeom>(topNode,TowerGeomNode + detector);
   if (!_epdgeom)
   {
      std::cout << "Could not locate geometry node " << TowerGeomNode + detector << std::endl;
      exit(1);
   }


   std::vector<EpHit> epdhitnorth;
   epdhitnorth.clear();

   std::vector<EpHit> epdhitsouth;
   epdhitsouth.clear();
      

   float tower_energy = 0.; float tile_phi = 0.;
   float eMax = 0.; float truncated_tile_e = 0.;
   int arm  = -1; int phibin = -1; int rbin = -1;
     
   unsigned int ntowers = _towerinfos->size();
   for (unsigned int ch = 0; ch < ntowers;  ch++)
   {

       TowerInfo *_tower = _towerinfos->get_tower_at_channel(ch);
       unsigned int key = TowerInfoDefs::encode_epd(ch);
       tower_energy = _tower->get_energy();

       if(tower_energy < 0.2) continue; //remove low response tiles
       
       tile_phi = _epdgeom->get_phi(key);
       arm = TowerInfoDefs::get_epd_arm(key);
       phibin = TowerInfoDefs::get_epd_phibin(key);
       rbin = TowerInfoDefs::get_epd_rbin(key);
       
       if(mTruncationInput)
       { 
         eMax = mTruncationInput->GetBinContent(arm+1,cent_index+1,rbin+1);
         std::cout<<arm+1<<"\t"<<cent_index+1<<"\t"<<rbin+1<<mTruncationInput->GetBinContent(arm+1,cent_index+1,rbin+1)<<"\t"<<eMax<<std::endl;
       }       
       else  
       {
         eMax = 100.;
       } //cant find truncation file, I will essentially not truncate

       truncated_tile_e = (tower_energy < eMax) ? tower_energy : eMax; //do tile energy truncation
      
       if(arm == 0) //south
       {

        EpHit newepdHit;
        newepdHit.nMip = truncated_tile_e;
        newepdHit.phi = tile_phi;
        newepdHit.iy = phibin;
        newepdHit.ix = rbin;
        newepdHit.wheel = arm;      
    
        if (rbin == 0)
        {
          newepdHit.sameRing = &Sepd_phi_list0[newepdHit.ix];
        }
        else
        {
          newepdHit.sameRing = &Sepd_phi_list[newepdHit.ix];
        }

       epdhitsouth.push_back(newepdHit); //store this ephit for south wheel
      }
       
      else if(arm == 1) //north
       {

         EpHit newepdHit;
         newepdHit.phi = tile_phi;
         newepdHit.nMip = truncated_tile_e;
         newepdHit.iy = phibin;
         newepdHit.ix = rbin;
         newepdHit.wheel = arm;

        if (rbin == 0)
        {
         newepdHit.sameRing = &Nepd_phi_list0[newepdHit.ix];
        }
        else
        {
          newepdHit.sameRing = &Nepd_phi_list[newepdHit.ix];
        }
      
       epdhitnorth.push_back(newepdHit);

      }
   }//end loop over sEPD Tower Info Container
    
 
   EpFinder_det[0]->ResultsEPD(epdhitsouth, 0, _EpInfo_det[0]);
   EpFinder_det[1]->ResultsEPD(epdhitnorth, 0, _EpInfo_det[1]);
    
   epdhitsouth.clear();
   epdhitnorth.clear();
   
 }
  
  return;

}

int EpFinderReco::GetNodes(PHCompositeNode *topNode)
{
 
  for (unsigned int i = 0; i < EventPlaneNodeName.size(); ++i)
  {
            
    _EpInfo_det[i] = findNode::getClass<EpInfo>(topNode, EventPlaneNodeName[i]);

    if (!_EpInfo_det[i])
    {
      std::cout << PHWHERE << ": Could not find node:"<< EventPlaneNodeName[i] << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
    
}





