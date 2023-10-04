#include "EventPlaneReco.h"

#include "Eventplaneinfo.h"
#include "EventplaneinfoMap.h"
#include "EventplaneinfoMapv1.h"
#include "Eventplaneinfov1.h"

#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>

#include <epd/EpdGeom.h>

#include <bbc/BbcPmtInfoContainerV1.h>
#include <bbc/BbcGeom.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <vector>  // for vector
#include <array>  // for array
#include <cfloat>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <set>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

EventPlaneReco::EventPlaneReco(const std::string &name)
  : SubsysReco(name)
, m_MaxOrder(3)
{
    south_rawpsi.resize(m_MaxOrder);
    south_qraw.resize(m_MaxOrder);
    north_rawpsi.resize(m_MaxOrder);
    north_qraw.resize(m_MaxOrder);
    
    for (auto &vec : south_qraw)
    {
        vec.resize(2);
    }
    
    for (auto &vec : north_qraw)
    {
        vec.resize(2);
    }
}

int EventPlaneReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "======================= EventPlaneReco::InitRun() =======================" << std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  return CreateNodes(topNode);
}

int EventPlaneReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "EventPlaneReco::process_event -- entered" << std::endl;
  }

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  EventplaneinfoMap *epmap = findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!epmap)
  {
    std::cout << PHWHERE << "::ERROR - cannot find EventplaneinfoMap" << std::endl;
    exit(-1);
  }

  TowerInfoContainer *epd_towerinfo = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EPD");
  if (!epd_towerinfo)
  {
    std::cout << PHWHERE << "::ERROR - cannot find TOWERINFO_CALIB_EPD" << std::endl;
    exit(-1);
  }

  EpdGeom *_epdgeom = findNode::getClass<EpdGeom>(topNode,"TOWERGEOM_EPD");
  if (!_epdgeom)
  {
    std::cout << PHWHERE << "::ERROR - cannot find TOWERGEOM_EPD" << std::endl;
    exit(-1);
  }
    
  /*
  //mbd container missing 
  BbcPmtInfoContainerV1 *m_bbcpmts = findNode::getClass<BbcPmtInfoContainerV1>(topNode, "BbcPmtInfoContainerV1");
  if (!m_bbcpmts)
  {
    std::std::cout << PHWHERE << ":: No BbcPmtInfoContainerV1!" << std::std::endl; exit(1);
  }
    
  if (m_bbcpmts)
  {
    for (int ipmt=0; ipmt<128; ipmt++)
    {
      BbcPmtInfoV1 *bbcpmt = m_bbcpmts->get_pmt( ipmt );
      int arm = ipmt/64;
      float mbd_e = bbcpmt->get_q();
      if(arm == 0)
      {
        std::std::cout<<"this is mbd south \t"<<mbd_e<<std::std::endl;
      }
      else if(arm == 1)
      {
        std::std::cout<<"this is mbd north \t"<<mbd_e<<std::std::endl;
      }
            
    }
 }
      
*/    
     if (epd_towerinfo)
     {
         if (Verbosity())
         {
             std::cout << "EventPlaneReco::process_event -  epd_towerinfo" << std::endl;
         }
         
         unsigned int ntowers = epd_towerinfo->size();
         for (unsigned int ch = 0; ch < ntowers;  ch++)
         {
             TowerInfo *_tower = epd_towerinfo->get_tower_at_channel(ch);
             unsigned int key = TowerInfoDefs::encode_epd(ch);
             float epd_e = _tower->get_energy();
             if(epd_e < 0.5) continue;
             float tile_phi = _epdgeom->get_phi(key);
             int arm = TowerInfoDefs::get_epd_arm(key);
             float truncated_e = (epd_e < 6.0)? epd_e:6.0;
             
             if(arm == 0)
             {
    
                for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
                 {
                     double Cosine = cos(tile_phi * (double) order);
                     double Sine = sin(tile_phi * (double) order);
                     south_qraw[order - 1][0] += truncated_e * Cosine; //south Qn,x
                     south_qraw[order - 1][1] += truncated_e * Sine; //south Qn,y
                }
             }
             else if(arm == 1)
             {
                
                 for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
                 {
                     double Cosine = cos(tile_phi * (double) order);
                     double Sine = sin(tile_phi * (double) order);
                     north_qraw[order - 1][0] += truncated_e * Cosine; //north Qn,x
                     north_qraw[order - 1][1] += truncated_e * Sine; //north Qn,y
                }
             }
         }
     }

    

    //get raw psi
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
        south_rawpsi[order - 1] = GetPsiInRange(south_qraw[order - 1][0], south_qraw[order - 1][1], order);
        north_rawpsi[order - 1] = GetPsiInRange(north_qraw[order - 1][0], north_qraw[order - 1][1], order);
    }
        
    
    ResetMe();

   if (epd_towerinfo)
   {
       Eventplaneinfo *ep = new Eventplaneinfov1(Eventplaneinfo::EPTYPE::EPDS);
       ep->set_id(epmap->size());
       
       for (unsigned int i = 0; i < south_rawpsi.size(); i++)
       {
           ep->set_psi_raw(i, south_rawpsi[i]);
       }
       
       ep->insert_ep_ids(Eventplaneinfo::EPDS, 1);
       epmap->insert(ep);
       
       if (Verbosity() > 1)
       {
          ep->identify();
       }
       
   }
    
    if (epd_towerinfo)
    {
        Eventplaneinfo *ep = new Eventplaneinfov1(Eventplaneinfo::EPTYPE::EPDN);
        ep->set_id(epmap->size());
        
        for (unsigned int i = 0; i < north_rawpsi.size(); i++)
        {
            ep->set_psi_raw(i, north_rawpsi[i]);
        }

        ep->insert_ep_ids(Eventplaneinfo::EPDN, 2);
        epmap->insert(ep);

        if (Verbosity() > 1)
        {
            ep->identify();
        }
   }

  if (Verbosity())
  {
      epmap->identify();
  }
    

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventPlaneReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the GLOBAL stuff under a sub-node directory
  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode)
  {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  // create the EventplaneinfoMap
  EventplaneinfoMap *eps = findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!eps)
  {
    eps = new EventplaneinfoMapv1();
    PHIODataNode<PHObject> *EpMapNode = new PHIODataNode<PHObject>(eps, "EventplaneinfoMap", "PHObject");
    globalNode->addNode(EpMapNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
    
}
double EventPlaneReco::GetPsiInRange(double Qx, double Qy, unsigned int order) const
{
  double temp;
  if ((Qx == 0.0) || (Qy == 0.0))
    temp = NAN;
  else
  {
    temp = atan2(Qy, Qx) / ((double) order);
    double AngleWrapAround = (2.0 * M_PI) / (double) order;
    if (temp < 0.0)
      temp += AngleWrapAround;
    else if (temp > AngleWrapAround)
      temp -= AngleWrapAround;
  }
  return temp;
    

}


void EventPlaneReco::ResetMe()
{
  for (auto &vec : south_qraw)
  {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &vec : north_qraw)
  {
    std::fill(vec.begin(), vec.end(), 0.);
  }
}
