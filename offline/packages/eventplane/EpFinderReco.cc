#include "EpFinderReco.h"

#include "EpFinder.h"
#include "EpInfo.h"    // for EpInfo
#include "EpInfov1.h"  // for EpInfo

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <epd/EpdGeomV1.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoDefs.h>

#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

//#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <filesystem>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>  // for TH2F
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
 
using namespace std;

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
    //for sEPD, provide truncation file if tower energies should be truncated
    if(!m_TruncationFileName.empty())
     {
       if (std::filesystem::exists(m_TruncationFileName))
       {
         std::ifstream truncate_Nmip;
         truncate_Nmip.open(m_TruncationFileName, std::ifstream::in);
         if (truncate_Nmip.is_open())
         {
           while (!truncate_Nmip.eof())
           {
               for (int k = 0; k < 2; k++)
                {
                   for (int i = 0; i < 10; i++)
                   {
                    for (int j = 0; j < 16; j++)
                    {
                        truncate_Nmip >> m_Epd_Trunc_e[k][i][j];
                        if (!std::isfinite(m_Epd_Trunc_e[k][i][j]))
                        {
                           std::cout << "Truncation value at rbin " << j
                           << ", centbin " << i << ", arm " << k << " in " << m_TruncationFileName
                           << " is not finite: " << m_Epd_Trunc_e[k][i][j] << std::endl;
                           gSystem->Exit(1);
                           exit(1);
                       }
                      else
                       {
                        if(Verbosity() >= 1) std::cout<<"Found Truncation File and will be using values \t"<< k << "\t" << i << "\t" << j << "\t" << m_Epd_Trunc_e[k][i][j] <<std::endl;
                      }
                    }
                  }
              }
           }
          truncate_Nmip.close();
         }
       }
     }
    return CreateNodes(topNode);
}

int EpFinderReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
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

   int cent_index = cent->get_centile(CentralityInfo::PROP::bimp)/10;
     
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
       
       if((m_Epd_Trunc_e[arm][cent_index][rbin]) == 0.)
       { 
         eMax = 100.;
       } //cant find truncation file, I will essentially not truncate
       else  
       {
         eMax = m_Epd_Trunc_e[arm][cent_index][rbin];
       }
       
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





