#include "caloTowerEmbed.h"

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoDefs.h>

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TFile.h>
#include <TProfile.h>

#include <algorithm>  // for min
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ostream
#include <stdexcept>  // for runtime_error

//____________________________________________________________________________..
caloTowerEmbed::caloTowerEmbed(const std::string &name)
  : SubsysReco(name)
{
  std::cout << "caloTowerEmbed::caloTowerEmbed(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
caloTowerEmbed::~caloTowerEmbed()
{
  std::cout << "caloTowerEmbed::~caloTowerEmbed() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int caloTowerEmbed::InitRun(PHCompositeNode *topNode)
{
  if (m_dettype == CaloTowerDefs::CEMC)
  {
    m_detector = "CEMC";
    m_templateName = "CEMC_TEMPLATE";
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
  {
    m_detector = "HCALIN";
    m_templateName = "IHCAL_TEMPLATE";
  }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
  {
    m_detector = "HCALOUT";
    m_templateName = "OHCAL_TEMPLATE";
  }
  else if (m_dettype == CaloTowerDefs::ZDC)
  {
    m_detector = "ZDC";
  }
  else if (m_dettype == CaloTowerDefs::SEPD)
  {
    m_detector = "SEPD";
  }

  try
  {
    CreateNodeTree(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_embedwaveform)
  {
    // the pulse shape the data's sub-threshold signal is given in ZS towers:
    // the reco template, peak moved to the post sample and normalized to 1 there
    std::string url = CDBInterface::instance()->getUrl(m_templateName);
    TFile *fin = TFile::Open(url.c_str());
    TProfile *h_template = (fin && fin->IsOpen()) ? dynamic_cast<TProfile *>(fin->Get("waveform_template")) : nullptr;
    if (!h_template)
    {
      std::cout << PHWHERE << " no waveform_template for " << m_templateName << " at " << url << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    const float peak = h_template->GetBinCenter(h_template->GetMaximumBin());
    const float norm = h_template->Interpolate(peak);
    m_template.assign(m_nsamples, 0);
    for (int j = 0; j < m_nsamples; j++)
    {
      m_template.at(j) = h_template->Interpolate(j - m_postsample + peak) / norm;
    }
    fin->Close();
    delete fin;

    // software ZS threshold of the fit of the embedded waveforms, as common/Calo_Fitting.C
    if (m_szs_threshold < 0)
    {
      if (m_dettype == CaloTowerDefs::CEMC)
      {
        m_szs_threshold = 60;
      }
      else if (m_dettype == CaloTowerDefs::HCALIN)
      {
        m_szs_threshold = 30;
      }
      else
      {
        m_szs_threshold = 15;
      }
    }
    m_neta = (m_dettype == CaloTowerDefs::CEMC) ? 96 : 24;
    m_bandsignal.assign(m_nsignalbands, 0);
    std::cout << Name() << ": software ZS threshold " << m_szs_threshold << " ADC, data signal estimated in "
              << m_nsignalbands << " |eta| bands" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTowerEmbed::signal_band(unsigned int channel) const
{
  const unsigned int key = (m_dettype == CaloTowerDefs::CEMC) ? TowerInfoDefs::encode_emcal(channel) : TowerInfoDefs::encode_hcal(channel);
  const int eta = static_cast<int>(TowerInfoDefs::getCaloTowerEtaBin(key));
  const int abseta = std::min(eta, m_neta - 1 - eta);  // 0 at the edges, m_neta / 2 - 1 at eta = 0
  return std::clamp(abseta * m_nsignalbands / (m_neta / 2), 0, m_nsignalbands - 1);
}

//____________________________________________________________________________..
void caloTowerEmbed::EstimateDataSignal()
{
  // mean sub-threshold data signal of the HW ZS towers of this event, per |eta| band:
  // <post - pre> of the data minus <X6 - X0> of the pedestal waveforms of the same towers,
  // so the pedestal shape cancels on average
  std::vector<double> sum(m_nsignalbands, 0);
  std::vector<int> n(m_nsignalbands, 0);
  for (unsigned int channel = 0; channel < _data_towers->size(); channel++)
  {
    TowerInfo *data = _data_towers->get_tower_at_channel(channel);
    TowerInfo *ped = m_PedestalContainer->get_tower_at_channel(channel);
    const float pre = data->get_waveform_value(0);
    const float post = data->get_waveform_value(1);
    if (!data->get_isZS() || data->get_isNotInstr() || pre == 0 || post == 0 || ped->get_waveform_value(0) == 0)
    {
      continue;
    }
    const int band = signal_band(channel);
    sum[band] += (post - pre) - (ped->get_waveform_value(m_postsample) - ped->get_waveform_value(0));
    ++n[band];
  }
  for (int b = 0; b < m_nsignalbands; b++)
  {
    m_bandsignal[b] = n[b] ? sum[b] / n[b] : 0;
  }
  if (Verbosity() > 1)
  {
    std::cout << Name() << ": mean data signal of the HW ZS towers per |eta| band (edge to center):";
    for (int b = 0; b < m_nsignalbands; b++)
    {
      std::cout << " " << m_bandsignal[b] << " (" << n[b] << ")";
    }
    std::cout << std::endl;
  }
}

//____________________________________________________________________________..
int caloTowerEmbed::process_event(PHCompositeNode * topNode)
{
  ++m_eventNumber;

  if (Verbosity())
  {
    std::cout << "event " << m_eventNumber << " working on " << m_detector << std::endl;
  }

  if (m_embedwaveform)
  {
    std::string ped_nodename = "PEDESTAL_" + m_detector;
    m_PedestalContainer = findNode::getClass<TowerInfoContainer>(topNode, ped_nodename);

    if (!m_PedestalContainer)
    {
      std::cout << PHWHERE << " " << ped_nodename << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    EstimateDataSignal();
  }

  unsigned int ntowers = _data_towers->size();
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {

    
    _sim_towers->get_tower_at_channel(channel)->set_status(_data_towers->get_tower_at_channel(channel)->get_status());

    TowerInfo *caloinfo_data = _data_towers->get_tower_at_channel(channel);
    TowerInfo *caloinfo_sim = _sim_towers->get_tower_at_channel(channel);

    // the embedded waveform is complete even where the data was zero suppressed: drop the
    // data's ZS flag, else the calibration applies the ZS cross calibration (post - pre to
    // template scale) to the template fit.  The fit sets it again for software ZS towers
    if (m_embedwaveform)
    {
      caloinfo_sim->set_isZS(false);
    }

    if (m_embedwaveform)
    {
      // when the data is not ZS-ed
      if (!(caloinfo_data->get_isZS() || caloinfo_data->get_isNotInstr()))
      {
        // here we really don't want the m_samples being greater than what data has!
        for (int j = 0; j < m_nsamples; j++)
        {
          // superpose the waveforms
          caloinfo_sim->set_waveform_value(j, caloinfo_data->get_waveform_value(j) + caloinfo_sim->get_waveform_value(j));
        }
      }
      else
      {
        // if this is ZS-ed or empty(like the packet is gone) the data only has pre (sample 0)
        // and post (sample 6). With t the reco template (t_6 = 1), ped a waveform from a
        // pedestal run and A = post - pre:
        //  - sim (x_6 - x_0) below the software ZS threshold of the fit:
        //      x_0 = pre, x_6 = post, x_j = pre + A t_j + (ped_j - <ped>)
        //    a tower that stays below the threshold keeps exactly the data's post - pre; one that
        //    goes above (data + sim) is fitted, with somewhat more than the real noise
        //  - sim alone at or above the threshold (the tower is fitted):
        //      x_j = pre + ped_j - ped_0 + s t_j
        //    the whole pedestal waveform gives the fit the real noise; the data's own post - pre
        //    is replaced by s, the mean data signal of the HW ZS towers of this event in this
        //    |eta| band (EstimateDataSignal)
        TowerInfo *pedestal_tower = m_PedestalContainer->get_tower_at_channel(channel);
        float pedestal_mean = 0;
        std::vector<float> m_waveform_pedestal;
        m_waveform_pedestal.resize(m_nsamples);
        for (int j = 0; j < m_nsamples; j++)
        {
          m_waveform_pedestal.at(j) = (j < m_datasamples) ? pedestal_tower->get_waveform_value(j) : pedestal_tower->get_waveform_value(m_datasamples - 1);
        }
        // <ped> over the samples the pedestal run has, not the padding
        const int nped = std::min(m_nsamples, m_datasamples);
        for (int j = 0; j < nped; j++)
        {
          pedestal_mean += m_waveform_pedestal.at(j);
        }
        pedestal_mean /= nped;

        float pre = caloinfo_data->get_waveform_value(0);
        float post = caloinfo_data->get_waveform_value(1);
        // no data (packet gone) or a lost post sample (the reco flags it with chi2 1e6):
        // take the baseline from the pedestal run and no data signal
        bool nodata = false;
        if (caloinfo_data->get_isNotInstr() || (pre != 0 && post == 0))
        {
          pre = pedestal_mean;
          post = pedestal_mean;
          nodata = true;
        }
        const float amplitude = post - pre;

        const float sim_amplitude = caloinfo_sim->get_waveform_value(m_postsample) - caloinfo_sim->get_waveform_value(0);
        if (sim_amplitude >= m_szs_threshold)
        {
          ++m_npedestalonly;
          const float s = nodata ? 0 : m_bandsignal[signal_band(channel)];
          for (int j = 0; j < m_nsamples; j++)
          {
            const float x = pre + m_waveform_pedestal.at(j) - m_waveform_pedestal.at(0) + s * m_template.at(j);
            caloinfo_sim->set_waveform_value(j, caloinfo_sim->get_waveform_value(j) + x);
          }
          continue;
        }

        for (int j = 0; j < m_nsamples; j++)
        {
          float x = pre + (amplitude * m_template.at(j)) + ((m_waveform_pedestal.at(j) - pedestal_mean) * m_pedestal_scale);
          if (j == 0)
          {
            x = pre;
          }
          else if (j == m_postsample)
          {
            x = post;
          }
          caloinfo_sim->set_waveform_value(j, caloinfo_sim->get_waveform_value(j) + x);
        }
      }
    }
    else
    {
      float data_E = caloinfo_data->get_energy();
      float sim_E = caloinfo_sim->get_energy();
      float embed_E = data_E + sim_E;

      _sim_towers->get_tower_at_channel(channel)->set_energy(embed_E);
      _sim_towers->get_tower_at_channel(channel)->set_time(caloinfo_data->get_time());
    }

  }  // end loop over channels

  return Fun4AllReturnCodes::EVENT_OK;
}

int caloTowerEmbed::End(PHCompositeNode * /*topNode*/)
{
  if (m_embedwaveform)
  {
    std::cout << Name() << ": " << m_npedestalonly << " HW ZS towers with sim above " << m_szs_threshold
              << " ADC took the pedestal waveform, " << m_npedestalonly / std::max(1., m_eventNumber + 1.) << " per event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void caloTowerEmbed::CreateNodeTree(PHCompositeNode *topNode)
{
  std::cout << "creating node" << std::endl;
  
  std::string TowerNodeName = m_inputNodePrefix + m_detector;
  std::string SimTowerNodeName = TowerNodeName;
  if (m_embedwaveform)
  {
    SimTowerNodeName = m_waveformNodePrefix + m_detector;
  }
  std::string GeomNodeName = "TOWERGEOM_" + m_detector;
  if (m_useRetower && m_detector == "CEMC")
  {
    TowerNodeName = m_inputNodePrefix + m_detector + "_RETOWER";
    GeomNodeName = "TOWERGEOM_HCALIN";
  }

  Fun4AllServer *se = Fun4AllServer::instance();
  PHCompositeNode *dataTopNode = se->topNode("TOPData");

  tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, GeomNodeName);
  if (!tower_geom)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "tower geom " << GeomNodeName << " missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find " + GeomNodeName + " node");
  }

  PHNodeIterator simIter(topNode);
  PHNodeIterator dataIter(dataTopNode);

  // data top node first

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(dataIter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerCalibration::CreateNodes");
  }

  PHCompositeNode *dstNodeSim = dynamic_cast<PHCompositeNode *>(simIter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNodeSim)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DSTsim Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DSTsim node in RawTowerCalibration::CreateNodes");
  }

  // data
  _data_towers = findNode::getClass<TowerInfoContainer>(dstNode, TowerNodeName);
  if (!_data_towers)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << TowerNodeName << " Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find " + TowerNodeName + " node in caloTowerEmbed::CreateNodes");
  }

  // sim
  _sim_towers = findNode::getClass<TowerInfoContainer>(dstNodeSim, SimTowerNodeName);
  if (!_sim_towers)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << TowerNodeName << " Sim Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find " + TowerNodeName + " Sim node in caloTowerEmbed::CreateNodes");
  }

  PHNodeIterator dstIterSim(dstNodeSim);
  PHCompositeNode *caloNode = dynamic_cast<PHCompositeNode *>(dstIterSim.findFirst("PHCompositeNode", m_detector));

  if (!caloNode)
  {
    caloNode = new PHCompositeNode(m_detector);
    dstNodeSim->addNode(caloNode);
  }

  return;
}
