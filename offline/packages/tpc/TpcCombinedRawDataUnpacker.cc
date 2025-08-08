#include "TpcCombinedRawDataUnpacker.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, hitsetkey
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <ffarawobjects/TpcRawHit.h>
#include <ffarawobjects/TpcRawHitContainer.h>

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdint>   // for exit
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, endl, bas...
#include <map>       // for _Rb_tree_iterator
#include <memory>
#include <utility>

TpcCombinedRawDataUnpacker::TpcCombinedRawDataUnpacker(std::string const& name, std::string const& outF)
  : SubsysReco(name)
  , outfile_name(outF)
{
  // Do nothing
}
TpcCombinedRawDataUnpacker::~TpcCombinedRawDataUnpacker()
{
  delete m_cdbttree;
}
void TpcCombinedRawDataUnpacker::ReadZeroSuppressedData()
{
  m_do_zs_emulation = true;
  m_do_baseline_corr = false;
  auto cdb = CDBInterface::instance();
  std::string dir = cdb->getUrl("TPC_ZS_THRESHOLDS");

  auto cdbtree = std::make_unique<CDBTTree>(dir);
  std::ostringstream name;
  if(dir.empty())
 {
  if(Verbosity() > 1)
  {
    std::cout << "use default tpc zs threshold of 20" << std::endl;
  }
  return;
 }
 
  for(int i=0; i<3; i++)
  {
    name.str("");
    name << "R"<<i+1<<"ADUThreshold";
    m_zs_threshold[i] = cdbtree->GetSingleFloatValue(name.str().c_str());
    if(Verbosity() > 1)
    {
      std::cout << "Loading ADU threshold of " << m_zs_threshold[i] << " for region " << i << std::endl;
    }
  }

}
int TpcCombinedRawDataUnpacker::Init(PHCompositeNode* /*topNode*/)
{
  std::cout << "TpcCombinedRawDataUnpacker::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  m_cdb = CDBInterface::instance();
  std::string calibdir = m_cdb->getUrl("TPC_FEE_CHANNEL_MAP");

  if (calibdir[0] == '/')
  {
    // use generic CDBTree to load
    m_cdbttree = new CDBTTree(calibdir);
    m_cdbttree->LoadCalibrations();
  }
  else
  {
    std::cout << "TpcRawDataDecoder::::InitRun No calibration file found" << std::endl;
    exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCombinedRawDataUnpacker::InitRun(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    std::cout << "TpcCombinedRawDataUnpacker::InitRun(PHCompositeNode* topNode)" << std::endl;
    std::cout << "\tCould not retrieve topNode; doing nothing" << std::endl;
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator dst_itr(topNode);
  PHCompositeNode* dst_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    if (Verbosity())
    {
      std::cout << "TpcCombinedRawDataUnpacker::InitRun(PHCompositeNode* topNode)" << std::endl;
    }
    if (Verbosity())
    {
      std::cout << "\tCould not retrieve dst_node; doing nothing" << std::endl;
    }
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator trkr_itr(dst_node);
  PHCompositeNode* trkr_node = dynamic_cast<PHCompositeNode*>(trkr_itr.findFirst("PHCompositeNode", "TRKR"));
  if (!trkr_node)
  {
    trkr_node = new PHCompositeNode("TRKR");
    dst_node->addNode(trkr_node);
  }

  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    if (Verbosity())
    {
      std::cout << "TpcCombinedRawDataUnpacker::InitRun(PHCompositeNode* topNode)" << std::endl;
    }
    if (Verbosity())
    {
      std::cout << "\tMaking TrkrHitSetContainer" << std::endl;
    }

    trkr_hit_set_container = new TrkrHitSetContainerv1;
    PHIODataNode<PHObject>* new_node = new PHIODataNode<PHObject>(trkr_hit_set_container, "TRKR_HITSET", "PHObject");
    trkr_node->addNode(new_node);
  }

  TpcRawHitContainer* tpccont = findNode::getClass<TpcRawHitContainer>(topNode, m_TpcRawNodeName);
  if (!tpccont)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "TpcCombinedRawDataUnpacker::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_TpcRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Removing module" << std::endl;

    Fun4AllServer* se = Fun4AllServer::instance();
    se->unregisterSubsystem(this);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (m_writeTree)
  {
    m_file = new TFile(outfile_name.c_str(), "RECREATE");
    m_ntup = new TNtuple("NT", "NT", "event:gtmbco:packid:ep:sector:side:fee:rx:entries:ped:width");
    m_ntup_hits = new TNtuple("NTH", "NTH", "event:gtmbco:packid:ep:sector:side:fee:chan:sampadd:sampch:phibin:tbin:layer:adc:ped:width");
    m_ntup_hits_corr = new TNtuple("NTC", "NTC", "event:gtmbco:packid:ep:sector:side:fee:chan:sampadd:sampch:phibin:tbin:layer:adc:ped:width:corr");
    if (m_doChanHitsCut)
    {
      m_HitChanDis = new TH2F("HitChanDis", "HitChanDis", 451, -0.5, 450.5, 256, -0.5, 255.5);
      m_HitsinChan = new TH1F("HitsinChan", "HitsinChan", 451, -0.5, 450.5);
    }
  }

  if (Verbosity() >= 1)
  {
    std::cout << "TpcCombinedRawDataUnpacker:: _ped_sig_cut = " << m_ped_sig_cut << std::endl;
    std::cout << "TpcCombinedRawDataUnpacker:: startevt = " << startevt << std::endl;
    std::cout << "TpcCombinedRawDataUnpacker:: endevt = " << endevt << std::endl;
  }

  // check run number if presamples need to be shifted, which went from 80 -> 120
  // at 41624
  Fun4AllServer* se = Fun4AllServer::instance();
  if (se->RunNumber() < 41624)
  {
    m_presampleShift = 0;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCombinedRawDataUnpacker::process_event(PHCompositeNode* topNode)
{
  if (_ievent < startevt || _ievent > endevt)
  {
    if (Verbosity() > 1)
    {
      std::cout << " Skip event " << _ievent << std::endl;
    }
    _ievent++;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  _ievent++;

  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "TpcCombinedRawDataUnpacker::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"TRKR_HITSET\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;
    gSystem->Exit(1);
    exit(1);

    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  TpcRawHitContainer* tpccont = findNode::getClass<TpcRawHitContainer>(topNode, m_TpcRawNodeName);
  if (!tpccont)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "TpcCombinedRawDataUnpacker::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_TpcRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;

    gSystem->Exit(1);
    exit(1);
  }

  PHG4TpcCylinderGeomContainer* geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  TrkrDefs::hitsetkey hit_set_key = 0;
  TrkrDefs::hitkey hit_key = 0;
  TrkrHitSetContainer::Iterator hit_set_container_itr;
  TrkrHit* hit = nullptr;

  uint64_t bco_min = UINT64_MAX;
  uint64_t bco_max = 0;

  const auto nhits = tpccont->get_nhits();

  int max_time_range = 0;

  for (unsigned int i = 0; i < nhits; i++)
  {
    TpcRawHit* tpchit = tpccont->get_hit(i);
    uint64_t gtm_bco = tpchit->get_gtm_bco();

    bco_min = std::min(gtm_bco, bco_min);
    bco_max = std::max(gtm_bco, bco_max);

    int fee = tpchit->get_fee();
    int channel = tpchit->get_channel();

    int feeM = FEE_map[fee];
    if (FEE_R[fee] == 2)
    {
      feeM += 6;
    }
    if (FEE_R[fee] == 3)
    {
      feeM += 14;
    }

    int side = 1;
    int32_t packet_id = tpchit->get_packetid();
    int ep = (packet_id - 4000) % 10;
    int sector = (packet_id - 4000 - ep) / 10;
    if (sector > 11)
    {
      side = 0;
    }

    unsigned int key = (256 * (feeM)) + channel;
    std::string varname = "layer";
    int layer = m_cdbttree->GetIntValue(key, varname);
    // antenna pads will be in 0 layer
    if (layer <= 6)
    {
      continue;
    }

    uint16_t sampadd = tpchit->get_sampaaddress();
    uint16_t sampch = tpchit->get_sampachannel();
    //    uint16_t sam = tpchit->get_samples();
    max_time_range = tpchit->get_samples();
    varname = "phi";  // + std::to_string(key);
    double phi = ((side == 1 ? 1 : -1) * (m_cdbttree->GetDoubleValue(key, varname) - M_PI / 2.)) + ((sector % 12) * M_PI / 6);
    PHG4TpcCylinderGeom* layergeom = geom_container->GetLayerCellGeom(layer);
    unsigned int phibin = layergeom->get_phibin(phi, side);
    unsigned int region = 0;
    if(layer > 15)
    {
      region = 1;
    }
    if( layer > 31)
    {
      region = 2;
    }
    hit_set_key = TpcDefs::genHitSetKey(layer, (mc_sectors[sector % 12]), side);
    hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);

    float hpedestal = 0;
    float hpedwidth = 0;

    if (Verbosity() > 2)
    {
      std::cout << "TpcCombinedRawDataUnpacker:: do zero suppression" << std::endl;
    }
    TH2* feehist = nullptr;
    hpedestal = 60;
    hpedwidth = m_zs_threshold[region];

    unsigned int pad_key = create_pad_key(side, layer, phibin);

    std::map<unsigned int, chan_info>::iterator chan_it = chan_map.find(pad_key);
    if (chan_it != chan_map.end())
    {
      (*chan_it).second.ped = hpedestal;
      (*chan_it).second.width = hpedwidth;
    }
    else
    {
      chan_info nucinfo;
      nucinfo.fee = fee;
      nucinfo.ped = hpedestal;
      nucinfo.width = hpedwidth;
      chan_map.insert(std::make_pair(pad_key, nucinfo));
    }
    int rx = get_rx(layer);
    unsigned int fee_key = create_fee_key(side, mc_sectors[sector % 12], rx, fee);
    // find or insert TH2C;
    std::map<unsigned int, TH2*>::iterator fee_map_it;

    fee_map_it = feeadc_map.find(fee_key);
    if (fee_map_it != feeadc_map.end())
    {
      feehist = (*fee_map_it).second;
    }
    else
    {
      std::string histname = "h" + std::to_string(fee_key);
      feehist = new TH2C(histname.c_str(), "histname", max_time_range + 1, -0.5, max_time_range + 0.5, 501, -0.5, 1000.5);
      feeadc_map.insert(std::make_pair(fee_key, feehist));
      std::vector<int> feeentries(feehist->GetNbinsX(), 0);
      feeentries_map.insert(std::make_pair(fee_key, feeentries));
    }
    auto fee_entries_it = feeentries_map.find(fee_key);
    std::vector<int>& fee_entries_vec = (*fee_entries_it).second;

    float threshold_cut = m_zs_threshold[region];

    int nhitschan = 0;

    if (m_doChanHitsCut)
    {
      for (std::unique_ptr<TpcRawHit::AdcIterator> adc_iterator(tpchit->CreateAdcIterator()); !adc_iterator->IsDone(); adc_iterator->Next())
      {
        const uint16_t s = adc_iterator->CurrentTimeBin();
        const uint16_t adc = adc_iterator->CurrentAdc();
        int t = s - m_presampleShift - m_t0;
        if (t < 0)
        {
          continue;
        }
        if (feehist != nullptr)
        {
          if (adc > 0)
          {
            if ((float(adc) - hpedestal) > threshold_cut)
            {
              nhitschan++;
            }
          }
        }
      }
      if (m_writeTree)
      {
        m_HitChanDis->Fill(nhitschan, channel);
        m_HitsinChan->Fill(nhitschan);
      }
      if (nhitschan > m_ChanHitsCut)
      {
        continue;
      }
    }

    for (std::unique_ptr<TpcRawHit::AdcIterator> adc_iterator(tpchit->CreateAdcIterator());
         !adc_iterator->IsDone();
         adc_iterator->Next())
    {
      const uint16_t s = adc_iterator->CurrentTimeBin();
      const uint16_t adc = adc_iterator->CurrentAdc();
      int t = s - m_presampleShift - m_t0;
      if (t < 0)
      {
        continue;
      }
      if (feehist != nullptr)
      {
        if (adc > 0)
        {
          if ((float(adc) - hpedestal) > threshold_cut)
          {
            feehist->Fill(t, adc - hpedestal);
            if (t < (int) fee_entries_vec.size())
            {
              fee_entries_vec[t]++;
            }
          }
        }
      }

      if ((float(adc) - hpedestal) > threshold_cut)
      {
        hit_key = TpcDefs::genHitKey(phibin, (unsigned int) t);
        // find existing hit, or create new one
        hit = hit_set_container_itr->second->getHit(hit_key);
        if (!hit)
        {
          hit = new TrkrHitv2();
          hit->setAdc(float(adc) - hpedestal);
          hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
        }

        if (m_writeTree)
        {
          float fXh[18];
          int nh = 0;

          fXh[nh++] = _ievent - 1;
          fXh[nh++] = gtm_bco;                  // gtm_bco;
          fXh[nh++] = packet_id;                // packet_id;
          fXh[nh++] = ep;                       // ep;
          fXh[nh++] = mc_sectors[sector % 12];  // Sector;
          fXh[nh++] = side;
          fXh[nh++] = fee;
          fXh[nh++] = channel;  // channel;
          fXh[nh++] = sampadd;  // sampadd;
          fXh[nh++] = sampch;   // sampch;
          fXh[nh++] = (float) phibin;
          fXh[nh++] = (float) t;
          fXh[nh++] = layer;
          fXh[nh++] = (float(adc) - hpedestal);
          fXh[nh++] = hpedestal;
          fXh[nh++] = hpedwidth;
          m_ntup_hits->Fill(fXh);
        }
      }
    }
  }

  if (m_do_baseline_corr == true)
  {
    // Histos filled now process them for fee local baselines

    int nhistfilled = 0;
    int nhisttotal = 0;
    for (auto& hiter : feeadc_map)
    {
      if (hiter.second != nullptr)
      {
        unsigned int fee_key = hiter.first;
        unsigned int side;
        unsigned int sector;
        unsigned int rx;
        unsigned int fee;
        unpack_fee_key(side, sector, rx, fee, fee_key);
        TH2* hist2d = hiter.second;
        std::map<unsigned int, std::vector<int>>::iterator fee_entries_it = feeentries_map.find(fee_key);
        if (fee_entries_it == feeentries_map.end())
        {
          continue;
          //	  fee_entries_vec_it = (*fee_entries_it).second.begin();
        }
        std::vector<int>::iterator fee_entries_vec_it = (*fee_entries_it).second.begin();

        std::vector<float> pedvec(hist2d->GetNbinsX(), 0);
        feebaseline_map.insert(std::make_pair(hiter.first, pedvec));
        std::map<unsigned int, std::vector<float>>::iterator fee_blm_it = feebaseline_map.find(hiter.first);
        (*fee_blm_it).second.resize(hist2d->GetNbinsX(), 0);
        for (int binx = 1; binx < hist2d->GetNbinsX(); binx++)
        {
          double timebin = (hist2d->GetXaxis())->GetBinCenter(binx);
          std::string histname1d = "h" + std::to_string(hiter.first) + "_" + std::to_string((int) timebin);
          nhisttotal++;
          float local_ped = 0;
          float local_width = 0;
          float entries = fee_entries_vec_it[timebin];
          if (fee_entries_vec_it[timebin] > 100)
          {
            nhistfilled++;

            TH1D* hist1d = hist2d->ProjectionY(histname1d.c_str(), binx, binx);
            if (hist1d->GetEntries() != fee_entries_vec_it[timebin])
            {
              std::cout << " vec " << fee_entries_vec_it[timebin]
                        << " hist " << hist1d->GetEntries()
                        << std::endl;
            }
            if (hist1d->GetEntries() > 10)
            {
              int maxbin = hist1d->GetMaximumBin();
              // calc peak position
              double hadc_sum = 0.0;
              double hibin_sum = 0.0;
              double hibin2_sum = 0.0;

              for (int isum = -3; isum <= 3; isum++)
              {
                float val = hist1d->GetBinContent(maxbin + isum);
                float center = hist1d->GetBinCenter(maxbin + isum);
                hibin_sum += center * val;
                hibin2_sum += center * center * val;
                hadc_sum += val;
              }
              local_ped = hibin_sum / hadc_sum;
              local_width = sqrt((hibin2_sum / hadc_sum) - (local_ped * local_ped));
            }
            delete hist1d;
          }
          (*fee_blm_it).second[(int) timebin] = local_ped + m_baseline_nsigma * local_width;

          if (m_writeTree)
          {
            float fXh[11];
            int nh = 0;

            fXh[nh++] = _ievent - 1;
            fXh[nh++] = 0;                        // gtm_bco;
            fXh[nh++] = 0;                        // packet_id;
            fXh[nh++] = 0;                        // ep;
            fXh[nh++] = mc_sectors[sector % 12];  // Sector;
            fXh[nh++] = side;
            fXh[nh++] = fee;
            fXh[nh++] = rx;
            fXh[nh++] = entries;
            fXh[nh++] = local_ped;
            fXh[nh++] = local_width;
            m_ntup->Fill(fXh);
          }
        }
      }
    }
    if (Verbosity() >= 1)
    {
      std::cout << " filled " << nhistfilled
                << " total " << nhisttotal
                << std::endl;

      std::cout << "second loop " << m_do_baseline_corr << std::endl;
    }

    // second loop over hits to apply baseline correction
    TrkrHitSetContainer::ConstRange hitsetrange;
    hitsetrange = trkr_hit_set_container->getHitSets(TrkrDefs::TrkrId::tpcId);

    for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
         hitsetitr != hitsetrange.second;
         ++hitsetitr)
    {
      // if(count>0)continue;
      TrkrHitSet* hitset = hitsetitr->second;
      unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
      int side = TpcDefs::getSide(hitsetitr->first);
      unsigned int sector = TpcDefs::getSectorId(hitsetitr->first);

      TrkrHitSet::ConstRange hitrangei = hitset->getHits();

      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
           hitr != hitrangei.second;
           ++hitr)
      {
        unsigned short phibin = TpcDefs::getPad(hitr->first);
        unsigned short tbin = TpcDefs::getTBin(hitr->first);
        unsigned short adc = (hitr->second->getAdc());

        unsigned int pad_key = create_pad_key(side, layer, phibin);

        float fee = 0;
        std::map<unsigned int, chan_info>::iterator chan_it = chan_map.find(pad_key);
        if (chan_it != chan_map.end())
        {
          chan_info cinfo = (*chan_it).second;
          fee = cinfo.fee;
          // hpedestal2 = cinfo.ped;
          // hpedwidth2 = cinfo.width;
        }

        int rx = get_rx(layer);
        float corr = 0;

        unsigned int fee_key = create_fee_key(side, sector, rx, fee);
        std::map<unsigned int, std::vector<float>>::iterator fee_blm_it = feebaseline_map.find(fee_key);
        if (fee_blm_it != feebaseline_map.end())
        {
          if (tbin < (int) (*fee_blm_it).second.size())
          {
            corr = (*fee_blm_it).second[tbin];
          }
          hitr->second->setAdc(0);
          float nuadc = (float(adc) - corr);
          nuadc = std::max<float>(nuadc, 0);
          hitr->second->setAdc(nuadc);

          if (m_writeTree)
          {
            float fXh[18];
            int nh = 0;

            fXh[nh++] = _ievent - 1;
            fXh[nh++] = 0;       // gtm_bco;
            fXh[nh++] = 0;       // packet_id;
            fXh[nh++] = 0;       // ep;
            fXh[nh++] = sector;  // mc_sectors[sector % 12];//Sector;
            fXh[nh++] = side;
            fXh[nh++] = fee;
            fXh[nh++] = 0;  // channel;
            fXh[nh++] = 0;  // sampadd;
            fXh[nh++] = 0;  // sampch;
            fXh[nh++] = (float) phibin;
            fXh[nh++] = (float) tbin;
            fXh[nh++] = layer;
            fXh[nh++] = float(adc);
            fXh[nh++] = 0;  // hpedestal2;
            fXh[nh++] = 0;  // hpedwidth2;
            fXh[nh++] = corr;

            m_ntup_hits_corr->Fill(fXh);
          }
        }
      }
    }
  }
  // reset histogramms
  for (auto& hiter2 : feeadc_map)
  {
    if (hiter2.second != nullptr)
    {
      hiter2.second->Reset();
    }
  }
  feebaseline_map.clear();
  for (auto& hiter_entries : feeentries_map)
  {
    hiter_entries.second.assign(hiter_entries.second.size(), 0);
  }

  if (Verbosity())
  {
    std::cout << " event BCO: " << bco_min << " - " << bco_max << std::endl;
    std::cout << "TpcCombinedRawDataUnpacker:: done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCombinedRawDataUnpacker::End(PHCompositeNode* /*topNode*/)
{
  if (m_writeTree)
  {
    m_file->cd();
    m_ntup->Write();
    m_ntup_hits->Write();
    m_ntup_hits_corr->Write();
    if (m_doChanHitsCut)
    {
      m_HitsinChan->Write();
    }
    m_file->Close();
  }
  if (Verbosity())
  {
    std::cout << "TpcCombinedRawDataUnpacker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  // if(m_Debug==1) hm->dumpHistos(m_filename, "RECREATE");

  return Fun4AllReturnCodes::EVENT_OK;
}
