
#include "CaloWaveformSim.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for hit_idbits
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/EventHeader.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerSimv2.h>

#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom_Spacalv1.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom_Spacalv1.h>  // for PHG4CylinderGeom_Spaca...
#include <g4detectors/PHG4CylinderGeom_Spacalv3.h>

#include <TF1.h>
#include <TFile.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TTree.h>
#include <cassert>
#include <sstream>
#include <string>

double CaloWaveformSim::template_function(double *x, double *par)
{
  Double_t v1 = par[0] * h_template->Interpolate(x[0] - par[1]) + par[2];
  return v1;
}

CaloWaveformSim::CaloWaveformSim(const std::string &name)
  : SubsysReco(name)
{
}

CaloWaveformSim::~CaloWaveformSim()
{
  gsl_rng_free(m_RandomGenerator);
}

int CaloWaveformSim::InitRun(PHCompositeNode *topNode)
{
  // Initialize random generator
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int seed = PHRandomSeed();  // fixed seed handled in PHRandomSeed()
  gsl_rng_set(m_RandomGenerator, seed);

  // Load waveform template
  const char *calibroot = getenv("CALIBRATIONROOT");
  if (!calibroot)
  {
    std::cerr << "CaloWaveformSim::InitRun missing CALIBRATIONROOT" << std::endl;
    exit(1);
  }
  std::string templatefilename = std::string(calibroot) + "/CaloWaveSim/" + m_templatefile;
  TFile *ft = TFile::Open(templatefilename.c_str());
  assert(ft && ft->IsOpen());
  h_template = static_cast<TProfile *>(ft->Get("hpwaveform"));

  // Determine run number
  EventHeader *evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  m_runNumber = evtHeader ? evtHeader->get_RunNumber() : -1;
  if (Verbosity() > 0)
  {
    std::cout << "CaloWaveformSim::InitRun Run Number: " << m_runNumber << std::endl;
  }

  // Detector-specific setup
  std::string url;
  if (m_dettype == CaloTowerDefs::CEMC)
  {
    m_detector = "CEMC";
    encode_tower = TowerInfoDefs::encode_emcal;
    decode_tower = TowerInfoDefs::decode_emcal;
    m_sampling_fraction = 2e-02;
    m_nchannels = 24576;
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
  {
    m_detector = "HCALIN";
    encode_tower = TowerInfoDefs::encode_hcal;
    decode_tower = TowerInfoDefs::decode_hcal;
    m_sampling_fraction = 0.162166;
    m_nchannels = 1536;
  }
  else  // HCALOUT
  {
    m_detector = "HCALOUT";
    encode_tower = TowerInfoDefs::encode_hcal;
    decode_tower = TowerInfoDefs::decode_hcal;
    m_sampling_fraction = 3.38021e-02;
    m_nchannels = 1536;
  }

  // Gain settings
  m_gain = m_highgain ? ((m_detector == "CEMC") ? 16 : 32) : 1;

  // Data energy calibration
  if (!m_overrideCalibName) m_calibName = m_detector + "_calib_ADC_to_ETower";
  if (!m_overrideFieldName) m_fieldname = m_detector + "_calib_ADC_to_ETower";
  url = m_giveDirectURL ? m_directURL : CDBInterface::instance()->getUrl(m_calibName);
  if (!url.empty())
    cdbttree = new CDBTTree(url);
  else
  {
    std::cerr << "CaloWaveformSim::InitRun No data calibration for " << m_calibName << std::endl;
    exit(1);
  }

  // MC energy calibration (optional)
  if (!m_overrideMCCalibName) m_MC_calibName = m_detector + "_MC_RECALIB";
  if (!m_overrideMCFieldName) m_MC_fieldname = m_detector + "_calib_ADC_to_ETower";
  url = m_giveDirectURL_MC ? m_directURL_MC : CDBInterface::instance()->getUrl(m_MC_calibName);
  if (!url.empty())
    cdbttree_MC = new CDBTTree(url);
  else if (Verbosity() > 0)
  {
    std::cout << "CaloWaveformSim::InitRun No MC calibration for " << m_MC_calibName << std::endl;
  }

  // Time calibration (data)
  if (!m_overrideTimeCalibName) m_calibName_time = m_detector + "_meanTime";
  if (m_giveDirectURL_time)
  {
    url = m_directURL_time;
  }
  else
  {
    url = CDBInterface::instance()->getUrl(m_calibName_time);
    if (url.empty())
    {
      if (m_dotimecalib)
      {
        std::cerr << "CaloWaveformSim::InitRun No time calibration for " << m_calibName_time << std::endl;
        exit(1);
      }
    }
  }
  if (m_dotimecalib) cdbttree_time = new CDBTTree(url);
  if (Verbosity() > 0 && m_dotimecalib)
    std::cout << "CaloWaveformSim::InitRun Time calibration from " << url << std::endl;

  // Time calibration (MC)
  if (!m_overrideMCTimeCalibName) m_MC_calibName_time = m_detector + "_MC_meanTime";
  if (m_giveDirectURL_MC_time)
  {
    url = m_directURL_MC_time;
    cdbttree_MC_time = new CDBTTree(url);
  }
  else
  {
    url = CDBInterface::instance()->getUrl(m_MC_calibName_time);
    if (!url.empty())
    {
      cdbttree_MC_time = new CDBTTree(url);
    }
    else if (m_dotimecalib)
    {
      std::cerr << "CaloWaveformSim::InitRun No MC time calibration for " << m_MC_calibName_time << std::endl;
      exit(1);
    }
  }

  // Prepare waveform buffers
  m_waveforms.assign(m_nchannels, std::vector<float>(m_nsamples));

  // Create node tree and finish
  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloWaveformSim::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "CaloWaveformSim::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  }
  // maybe we really need to get the geometry node in in every event(otherwise layergeom become invalid when we get to the second file in the list?):
  if (m_dettype == CaloTowerDefs::CEMC)
  {
    PHG4CylinderGeomContainer *layergeo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_CEMC");
    if (!layergeo)
    {
      std::cout << PHWHERE << " CYLINDERGEOM_CEMC Node missing, doing nothing." << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
    const PHG4CylinderGeom *layergeom_raw = layergeo->GetFirstLayerGeom();
    assert(layergeom_raw);

    layergeom = dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);
    assert(layergeom);

    PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_CEMC");
    if (!seggeo)
    {
      std::cout << PHWHERE << " CYLINDERCELLGEOM_CEMC Node missing, doing nothing." << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
    PHG4CylinderCellGeom *geo_raw = seggeo->GetFirstLayerCellGeom();
    geo = dynamic_cast<PHG4CylinderCellGeom_Spacalv1 *>(geo_raw);
  }

  // initialize the waveform
  for (auto &waveform : m_waveforms)
  {
    for (auto &sample : waveform)
    {
      sample = 0.;
    }
  }
  // waveform TH1
  TF1 *f_fit = new TF1(
      "f_fit", [this](double *x, double *par)
      { return this->template_function(x, par); },
      0, m_nsamples, 3);
  f_fit->SetParameter(0, 1.0);
  f_fit->SetParameter(1, 0.0);
  float template_peak = f_fit->GetMaximumX();
  float shift_of_shift = m_timeshiftwidth * gsl_rng_uniform(m_RandomGenerator);

  float _shiftval = m_peakpos + shift_of_shift - template_peak;

  f_fit->SetParameters(1, _shiftval, 0);

  // get G4Hits
  std::string nodename = "G4HIT_" + m_detector;
  PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
  if (!hits)
  {
    std::cout << PHWHERE << " " << nodename << " Node missing, doing nothing." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  // loop over hits
  for (PHG4HitContainer::ConstIterator hititer = hits->getHits().first; hititer != hits->getHits().second; hititer++)
  {
    PHG4Hit *hit = hititer->second;
    if (hit->get_t(1) - hit->get_t(0) > m_deltaT)
    {
      continue;
    }

    // timing cut

    // get eta phi bin
    unsigned short etabin = 0;
    unsigned short phibin = 0;
    float correction = 1.;
    maphitetaphi(hit, etabin, phibin, correction);
    unsigned int key = encode_tower(etabin, phibin);
    float calibconst = cdbttree->GetFloatValue(key, m_fieldname);
    float e_vis = hit->get_light_yield();
    e_vis *= correction;
    float e_dep = e_vis / m_sampling_fraction;
    float ADC = (calibconst != 0) ? e_dep / calibconst : 0.;
    ADC *= m_gain;

    if(cdbttree_MC)
    {
      float MC_calibconst = cdbttree_MC->GetFloatValue(key, m_MC_fieldname);
      ADC *= MC_calibconst;
    }

    // if we have a tower by tower mean time, we shift the simulated waveform peak to that accordingly
    if (m_dotimecalib)
    { 
      float meantime = cdbttree_time->GetFloatValue(key, m_fieldname_time);
      float MCmeantime = cdbttree_MC_time->GetFloatValue(key, m_MC_fieldname_time);
      assert(m_peakpos == 6); // the MC mean time is derived when m_peakpos is set to 6
      _shiftval = m_peakpos + shift_of_shift - template_peak + meantime - MCmeantime;
    }

    float t0 = hit->get_t(0) / m_sampletime;
    unsigned int tower_index = decode_tower(key);
    // here I will add the truth matching part
    //  for the cell reco, the truth matching info relys on edep not light yield, I will be consistent here :)
    TowerInfo *tower = m_CaloWaveformContainer->get_tower_at_channel(tower_index);
    TowerInfo::EdepMap &edepMap = tower->get_hitEdepMap();
    TowerInfo::ShowerEdepMap &showerMap = tower->get_showerEdepMap();

    int showerID = hit->get_shower_id();
    float hitEdep = hit->get_edep();

    edepMap[hit->get_hit_id()] += hitEdep;
    showerMap[showerID] += hitEdep;

    f_fit->SetParameters(ADC, _shiftval + t0, 0.);
    for (int i = 0; i < m_nsamples; i++)
    {
      m_waveforms.at(tower_index).at(i) += f_fit->Eval(i);
    }
  }

  // do noise here and add to waveform

  if (m_noiseType == NoiseType::NOISE_TREE)
  {
    std::string ped_nodename = "PEDESTAL_" + m_detector;
    m_PedestalContainer = findNode::getClass<TowerInfoContainer>(topNode, ped_nodename);

    if (!m_PedestalContainer)
    {
      std::cout << PHWHERE << " " << ped_nodename << " Node missing, doing nothing." << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }

  for (int i = 0; i < m_nchannels; i++)
  {
    std::vector<float> m_waveform_pedestal;
    m_waveform_pedestal.resize(m_nsamples);
    if (m_noiseType == NoiseType::NOISE_TREE)
    {
      TowerInfo *pedestal_tower = m_PedestalContainer->get_tower_at_channel(i);
      float pedestal_mean = 0;
      for (int j = 0; j < m_nsamples; j++)
      {
        m_waveform_pedestal.at(j) = (j < m_pedestalsamples) ? pedestal_tower->get_waveform_value(j) : pedestal_tower->get_waveform_value(m_pedestalsamples - 1);
        pedestal_mean += m_waveform_pedestal.at(j);
      }
      pedestal_mean /= m_nsamples;
      for (int j = 0; j < m_nsamples; j++)
      {
        m_waveform_pedestal.at(j) = (m_waveform_pedestal.at(j) - pedestal_mean) * m_pedestal_scale + pedestal_mean;
      }
    }
    for (int j = 0; j < m_nsamples; j++)
    {
      if (m_noiseType == NoiseType::NOISE_TREE)
      {
        // TowerInfo *pedestal_tower = m_PedestalContainer->get_tower_at_channel(i);
        // m_waveforms.at(i).at(j) += (j < m_pedestalsamples) ? pedestal_tower->get_waveform_value(j) : pedestal_tower->get_waveform_value(m_pedestalsamples - 1);
        m_waveforms.at(i).at(j) += m_waveform_pedestal.at(j);
      }
      if (m_noiseType == NoiseType::NOISE_GAUSSIAN)
      {
        m_waveforms.at(i).at(j) += gsl_ran_gaussian(m_RandomGenerator, m_gaussian_noise);
      }
      if (m_noiseType == NoiseType::NOISE_NONE)
      {
        m_waveforms.at(i).at(j) += m_fixpedestal;
      }
      // saturate at 2^14 - 1
      if (m_waveforms.at(i).at(j) > 16383)
      {
        m_waveforms.at(i).at(j) = 16383;
      }
      if (m_waveforms.at(i).at(j) < 0)
      {
        m_waveforms.at(i).at(j) = 0;
      }

      m_CaloWaveformContainer->get_tower_at_channel(i)->set_waveform_value(j, m_waveforms.at(i).at(j));
    }
  }
  delete f_fit;
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloWaveformSim::maphitetaphi(PHG4Hit *g4hit, unsigned short &etabin, unsigned short &phibin, float &correction)
{
  if (m_dettype == CaloTowerDefs::CEMC)
  {
    int scint_id = g4hit->get_scint_id();
    PHG4CylinderGeom_Spacalv3::scint_id_coder decoder(scint_id);
    std::pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(decoder.tower_ID, decoder.sector_ID);
    const int &tower_ID_z = tower_z_phi_ID.first;
    const int &tower_ID_phi = tower_z_phi_ID.second;
    PHG4CylinderGeom_Spacalv3::tower_map_t::const_iterator it_tower =
        layergeom->get_sector_tower_map().find(decoder.tower_ID);
    assert(it_tower != layergeom->get_sector_tower_map().end());

    const int etabin_cell = geo->get_etabin_block(tower_ID_z);  // block eta bin
    const int sub_tower_ID_x = it_tower->second.get_sub_tower_ID_x(decoder.fiber_ID);
    const int sub_tower_ID_y = it_tower->second.get_sub_tower_ID_y(decoder.fiber_ID);
    unsigned short etabinshort = etabin_cell * layergeom->get_n_subtower_eta() + sub_tower_ID_y;
    unsigned short phibin_cell = tower_ID_phi * layergeom->get_n_subtower_phi() + sub_tower_ID_x;
    etabin = etabinshort;
    phibin = phibin_cell;

    // correction for emcal fiber
    if (light_collection_model.use_fiber_model())
    {
      const double z = 0.5 * (g4hit->get_local_z(0) + g4hit->get_local_z(1));
      assert(not std::isnan(z));
      correction *= light_collection_model.get_fiber_transmission(z);
    }

    if (light_collection_model.use_fiber_model())
    {
      const double x = it_tower->second.get_position_fraction_x_in_sub_tower(decoder.fiber_ID);
      const double y = it_tower->second.get_position_fraction_y_in_sub_tower(decoder.fiber_ID);
      correction *= light_collection_model.get_light_guide_efficiency(x, y);
    }
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
  {
    // int layer = (g4hit->get_hit_id() >> PHG4HitDefs::hit_idbits);
    unsigned int iphi = (unsigned int) (g4hit->get_hit_id() >> PHG4HitDefs::hit_idbits) / 4;
    unsigned int ieta = g4hit->get_scint_id();

    etabin = ieta;
    phibin = iphi;
  }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
  {
    // int layer = (g4hit->get_hit_id() >> PHG4HitDefs::hit_idbits);
    unsigned int iphi = (unsigned int) (g4hit->get_hit_id() >> PHG4HitDefs::hit_idbits) / 5;
    unsigned int ieta = g4hit->get_scint_id();

    etabin = ieta;
    phibin = iphi;
  }
  else
  {
    std::cout << PHWHERE << " Invalid detector type " << m_dettype << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
}

//____________________________________________________________________________..
int CaloWaveformSim::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "CaloWaveformSim::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloWaveformSim::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator topNodeItr(topNode);
  // DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(topNodeItr.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "PHComposite node created: DST" << std::endl;
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator nodeItr(dstNode);
  PHCompositeNode *DetNode;
  // enum CaloTowerDefs::DetectorSystem and TowerInfoContainer::DETECTOR are different!!!!
  TowerInfoContainer::DETECTOR DetectorEnum = TowerInfoContainer::DETECTOR::DETECTOR_INVALID;
  std::string DetectorNodeName;

  if (m_dettype == CaloTowerDefs::CEMC)
  {
    DetectorEnum = TowerInfoContainer::DETECTOR::EMCAL;
    DetectorNodeName = "CEMC";
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
  {
    DetectorEnum = TowerInfoContainer::DETECTOR::HCAL;
    DetectorNodeName = "HCALIN";
  }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
  {
    DetectorEnum = TowerInfoContainer::DETECTOR::HCAL;
    DetectorNodeName = "HCALOUT";
  }
  else
  {
    std::cout << PHWHERE << " Invalid detector type " << m_dettype << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  DetNode = dynamic_cast<PHCompositeNode *>(nodeItr.findFirst("PHCompositeNode", DetectorNodeName));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(DetectorNodeName);
    dstNode->addNode(DetNode);
  }
  m_CaloWaveformContainer = new TowerInfoContainerSimv2(DetectorEnum);

  PHIODataNode<PHObject> *newTowerNode = new PHIODataNode<PHObject>(m_CaloWaveformContainer, "WAVEFORM_" + m_detector, "PHObject");
  DetNode->addNode(newTowerNode);
}
