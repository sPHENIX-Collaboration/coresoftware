#include "ZdcReco.h"
#include "Zdcinfo.h"
#include "Zdcinfov1.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <cdbobjects/CDBTTree.h>  // for CDBTTree
#include <ffamodules/CDBInterface.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <array>  // for array
#include <cfloat>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <set>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>   // for vector

ZdcReco::ZdcReco(const std::string &name)
  : SubsysReco(name)
{
}

int ZdcReco::InitRun(PHCompositeNode *topNode)
{
    if (!m_overrideCalibName)
    {
      m_calibName = "data_driven_zdc_calib";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "zdc_calib";
    }
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "ZdcReco::::InitRun No calibration file for domain " << m_calibName << " found" << std::endl;
      exit(1);
    }

    
    PHNodeIterator node_itr(topNode);
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(node_itr.findFirst("PHCompositeNode", "RUN"));
    if (!runNode)
    {
       std::cout << PHWHERE << "RUN Node not found - that is fatal" << std::endl;
       gSystem->Exit(1);
       exit(1);
    }

    PHNodeIterator runiter(runNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(runiter.findFirst("PHCompositeNode", m_Detector));
    if (!DetNode)
    {
       DetNode = new PHCompositeNode(m_Detector);
       runNode->addNode(DetNode);
    }

    m_zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
    if (!m_zdcinfo)
    {
       m_zdcinfo = new Zdcinfov1();
       PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_zdcinfo, "Zdcinfo", "PHObject");
       DetNode->addNode(newNode);
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int ZdcReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "ZdcReco::process_event -- entered" << std::endl;
  }

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

    ResetMe(); 
    radius_north = 0.;
    radius_south = 0.;
    _sumS = 0.;
    _sumN = 0.;
 
    TowerInfoContainer *zdc_towerinfo = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
    if (!zdc_towerinfo)
    {
      std::cout << PHWHERE << "::ERROR - cannot find TOWERS_ZDC" << std::endl;
      exit(-1);
    }
   
    if (zdc_towerinfo)
    {
        if (Verbosity())
        {
            std::cout << "ZdcReco::process_event -  zdc_towerinfo" << std::endl;
        }
        
        unsigned int ntowers = zdc_towerinfo->size();
        
        if(ntowers != 52)
        {
            std::cout << "ZdcReco::process_event -  zdc size mismatch" << std::endl;
        }
                
        //get smd info
        for (unsigned int ch = 0; ch < ntowers; ch++)
        {
            TowerInfo *_tower = zdc_towerinfo->get_tower_at_channel(ch);
            float zdc_e = _tower->get_energy();
            float zdc_time = _tower->get_time_float();
            
            if(TowerInfoDefs::isSMD(ch))
            {
                vsmdadc.push_back(zdc_e);
                vsmdtime.push_back(zdc_time);
            }
        }
        
        //check smd mapping
        int ssize = vsmdadc.size();
        if(ssize != 32)
        {
            std::cout<< "smd channel mapping error" << std::endl;
            if(vsmdtime.size() != 32)
            exit(1);
        }
        
        //apply time cuts per smd side
        for (int j = 0; j < ssize; j++)
        {
            if(j < 16)
            {
                if(vsmdtime[j] > 9.0 && vsmdtime[j] < 14.0)
                {
                    smd_adc[j] = vsmdadc[j];
                }
                else
                {
                    smd_adc[j] = 0.0;
                }
            }
            else
            {
                if(vsmdtime[j] > 6.0 && vsmdtime[j] < 12.0)
                {
                    smd_adc[j] = vsmdadc[j];
                }
                else
                {
                    smd_adc[j] = 0.0;
                }
            }
        }
        
        //get smd position
        CompSmdPos();
        
        radius_north = sqrt(smd_pos[1]*smd_pos[1] + smd_pos[0]*smd_pos[0]);
        radius_south = sqrt(smd_pos[3]*smd_pos[3] + smd_pos[2]*smd_pos[2]);
        
        //get zdc info
        for (unsigned int ch = 0; ch < ntowers; ch++)
        {
            TowerInfo *_tower = zdc_towerinfo->get_tower_at_channel(ch);
            float zdc_e = _tower->get_energy();
            float zdc_time = _tower->get_time_float();
            
            if(TowerInfoDefs::isZDC(ch))
            {
                vzdcadc.push_back(zdc_e);
                vzdctime.push_back(zdc_time);
            }
        }

        //check zdc mapping
        int zsize = vzdcadc.size();
        if(zsize != 16)
        {
            std::cout<< "zdc channel mapping error" << std::endl;
            if(vzdctime.size() != 16)
            exit(1);
        }
             
        //apply time cuts per zdc and get sums
        for (int i = 0; i < zsize; i++)
        {
            unsigned int key = TowerInfoDefs::encode_zdc(i);
            int arm = TowerInfoDefs::get_zdc_side(key);
            
            if(vzdctime[i] > 5.0 && vzdctime[i] < 9.0)
            {
                if(arm == 0)
                {
                    if(vzdcadc[0] > 65.0 && vzdcadc[2] > 20.0 && radius_south < 2.0)
                    {
                        _sumS = vzdcadc[0]*cdbttree->GetFloatValue(0, m_fieldname) + vzdcadc[2]*cdbttree->GetFloatValue(2, m_fieldname) + vzdcadc[4]*cdbttree->GetFloatValue(4, m_fieldname);
                    }
                }
                else if(arm == 1)
                {
                    if(vzdcadc[8] > 65.0 && vzdcadc[10] > 20.0 && radius_north < 2.0)
                    {
                        _sumN = vzdcadc[8]*cdbttree->GetFloatValue(8, m_fieldname) + vzdcadc[10]*cdbttree->GetFloatValue(10, m_fieldname) + vzdcadc[12]*cdbttree->GetFloatValue(12, m_fieldname);
                    }
                }
                
            }
        }
       
        m_zdcinfo->set_zdc_energy(0, _sumS);
        m_zdcinfo->set_zdc_energy(1, _sumN);
        
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

void ZdcReco::ResetMe()
{
  vsmdadc.clear();
  vsmdtime.clear();
  vzdcadc.clear();
  vzdctime.clear();
}

int ZdcReco::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
    
void ZdcReco::CompSmdPos()  //computing position with weighted averages
{
    float w_ave[4];  // 0 -> north hor; 1 -> noth vert; 2 -> south hor; 3 -> south vert.
    float weights[4] = {0};
    memset(weights, 0, sizeof(weights));  // memset float works only for 0
    float w_sum[4];
    memset(w_sum, 0, sizeof(w_sum));

    // these constants convert the SMD channel number into real dimensions (cm's)
    const float hor_scale = 2.0 * 11.0 / 10.5 * sin(M_PI / 4);  // from gsl_math.h
    const float ver_scale = 1.5 * 11.0 / 10.5;
    float hor_offset = (hor_scale * 8 / 2.0) * (7.0 / 8.0);
    float ver_offset = (ver_scale * 7 / 2.0) * (6.0 / 7.0);

    for (int i = 0; i < 8; i++)
    {
      weights[0] += smd_adc[i];  //summing weights
      weights[2] += smd_adc[i + 16];
      w_sum[0] += (float) i * smd_adc[i];  //summing for the average
      w_sum[2] += ((float) i + 16.) * smd_adc[i + 16];
    }
    for (int i = 0; i < 7; i++)
    {
      weights[1] += smd_adc[i + 8];
      weights[3] += smd_adc[i + 24];
      w_sum[1] += ((float) i + 8.) * smd_adc[i + 8];
      w_sum[3] += ((float) i + 24.) * smd_adc[i + 24];
    }
      
    if (weights[0] > 0.0)
    {
      w_ave[0] = w_sum[0] / weights[0];  //average = sum / sumn of weights...
      smd_pos[0] = hor_scale * w_ave[0] - hor_offset;
    }
    else
    {
      smd_pos[0] = 0;
    }
    if (weights[1] > 0.0)
    {
      w_ave[1] = w_sum[1] / weights[1];
      smd_pos[1] = ver_scale * (w_ave[1] - 8.0) - ver_offset;
    }
    else
    {
     smd_pos[1] = 0;
    }

    if (weights[2] > 0.0)
    {
      w_ave[2] = w_sum[2] / weights[2];
      smd_pos[2] = hor_scale * (w_ave[2] - 16.0) - hor_offset;
    }
    else
    {
      smd_pos[2] = 0;
    }

    if (weights[3] > 0.0)
    {
      w_ave[3] = w_sum[3] / weights[3];
      smd_pos[3] = ver_scale * (w_ave[3] - 24.0) - ver_offset;
    }
    else
    {
      smd_pos[3] = 0;
    }
}
