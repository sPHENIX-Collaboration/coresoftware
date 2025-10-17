#include "TpcCombinedRawDataUnpackerDebug.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, hitsetkey
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

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

#include <cmath>
#include <cstdint>   // for exit
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, endl, bas...
#include <map>       // for _Rb_tree_iterator
#include <memory>
#include <utility>

#define dEBUG

TpcCombinedRawDataUnpackerDebug::TpcCombinedRawDataUnpackerDebug(std::string const& name, std::string const& outF)
  : SubsysReco(name)
  , outfile_name(outF)
{
  // Do nothing
}

int TpcCombinedRawDataUnpackerDebug::Init(PHCompositeNode* /*topNode*/)
{
  std::cout << "TpcCombinedRawDataUnpackerDebug::Init(PHCompositeNode *topNode) Initializing" << std::endl;

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

int TpcCombinedRawDataUnpackerDebug::InitRun(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    std::cout << "TpcCombinedRawDataUnpackerDebug::InitRun(PHCompositeNode* topNode)" << std::endl;
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
      std::cout << "TpcCombinedRawDataUnpackerDebug::InitRun(PHCompositeNode* topNode)" << std::endl;
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
      std::cout << "TpcCombinedRawDataUnpackerDebug::InitRun(PHCompositeNode* topNode)" << std::endl;
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
    std::cout << "TpcCombinedRawDataUnpackerDebug::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_TpcRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Removing module" << std::endl;

    Fun4AllServer* se = Fun4AllServer::instance();
    se->unregisterSubsystem(this);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (m_writeTree)
  {
    m_file = new TFile(outfile_name.c_str(), "RECREATE");
    m_ntup = new TNtuple("NT", "NT", "event:gtmbco:packid:ep:sector:side:fee:chan:sampadd:sampch:nsamples");
    m_ntup_hits = new TNtuple("NTH", "NTH", "event:gtmbco:packid:ep:sector:side:fee:chan:sampadd:sampch:phibin:tbin:layer:adc:ped:width");
    m_ntup_hits_corr = new TNtuple("NTC", "NTC", "event:gtmbco:packid:ep:sector:side:fee:chan:sampadd:sampch:phibin:tbin:layer:adc:ped:width:corr");
  }

  if (Verbosity() >= 1)
  {
    std::cout << "TpcCombinedRawDataUnpackerDebug:: _do_zerosup = " << m_do_zerosup << std::endl;
    std::cout << "TpcCombinedRawDataUnpackerDebug:: _do_noise_rejection = " << m_do_noise_rejection << std::endl;
    std::cout << "TpcCombinedRawDataUnpackerDebug:: _ped_sig_cut = " << m_ped_sig_cut << std::endl;
    std::cout << "TpcCombinedRawDataUnpackerDebug:: startevt = " << startevt << std::endl;
    std::cout << "TpcCombinedRawDataUnpackerDebug:: endevt = " << endevt << std::endl;
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

int TpcCombinedRawDataUnpackerDebug::process_event(PHCompositeNode* topNode)
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
  TH1F pedhist("pedhist", "pedhist", 251, -2.0, 1002);

  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "TpcCombinedRawDataUnpackerDebug::process_event(PHCompositeNode* topNode)" << std::endl;
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
    std::cout << "TpcCombinedRawDataUnpackerDebug::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_TpcRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;

    gSystem->Exit(1);
    exit(1);
  }

  PHG4TpcGeomContainer* geom_container =
      findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TPCGEOMCONTAINER" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  TrkrDefs::hitsetkey hit_set_key = 0;
  TrkrDefs::hitkey hit_key = 0;
  TrkrHitSetContainer::Iterator hit_set_container_itr;
  TrkrHit* hit = nullptr;

  uint64_t bco_min = UINT64_MAX;
  uint64_t bco_max = 0;

  const auto nhits = tpccont->get_nhits();

  int ntotalchannels = 0;
  int n_noisychannels = 0;
  int max_time_range = 0;
  for (unsigned int i = 0; i < nhits; i++)
  {
    TpcRawHit* tpchit = tpccont->get_hit(i);
    uint64_t gtm_bco = tpchit->get_gtm_bco();

    if (gtm_bco < bco_min)
    {
      bco_min = gtm_bco;
    }
    if (gtm_bco > bco_max)
    {
      bco_max = gtm_bco;
    }

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

    unsigned int key = 256 * (feeM) + channel;
    std::string varname = "layer";
    int layer = m_cdbttree->GetIntValue(key, varname);
    // antenna pads will be in 0 layer
    if (layer <= 0)
    {
      continue;
    }

    uint16_t sampadd = tpchit->get_sampaaddress();
    uint16_t sampch = tpchit->get_sampachannel();
    uint16_t sam = tpchit->get_samples();
    max_time_range = sam;
    varname = "phi";  // + std::to_string(key);
    double phi = -1 * pow(-1, side) * m_cdbttree->GetDoubleValue(key, varname) - M_PI/2. + (sector % 12) * M_PI / 6;
    PHG4TpcGeom* layergeom = geom_container->GetLayerCellGeom(layer);
    unsigned int phibin = layergeom->get_phibin(phi, side);
    if (m_writeTree)
    {
      float fX[12];
      int n = 0;

      fX[n++] = _ievent - 1;
      fX[n++] = gtm_bco;
      fX[n++] = packet_id;
      fX[n++] = ep;
      fX[n++] = sector;
      fX[n++] = side;
      fX[n++] = fee;
      fX[n++] = channel;
      fX[n++] = sampadd;
      fX[n++] = sampch;
      fX[n++] = sam;
      m_ntup->Fill(fX);
    }

    hit_set_key = TpcDefs::genHitSetKey(layer, (mc_sectors[sector % 12]), side);
    hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);

    float hpedestal = 0;
    float hpedwidth = 0;
    pedhist.Reset();

    if (!m_do_zerosup)
    {
      if (Verbosity() > 2)
      {
        std::cout << "TpcCombinedRawDataUnpackerDebug:: no zero suppression" << std::endl;
      }
      for (std::unique_ptr<TpcRawHit::AdcIterator> adc_iterator(tpchit->CreateAdcIterator());
           !adc_iterator->IsDone();
           adc_iterator->Next())
      {
        const uint16_t s = adc_iterator->CurrentTimeBin();
        const uint16_t adc = adc_iterator->CurrentAdc();

        int t = s - m_presampleShift;

        hit_key = TpcDefs::genHitKey(phibin, (unsigned int) t);
        // find existing hit, or create new one
        hit = hit_set_container_itr->second->getHit(hit_key);
        if (!hit)
        {
          hit = new TrkrHitv2();
          hit->setAdc(float(adc));

          hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
        }
      }
    }
    else
    {
      if (Verbosity() > 2)
      {
        std::cout << "TpcCombinedRawDataUnpackerDebug:: do zero suppression" << std::endl;
      }
      TH2I* feehist = nullptr;
      if (!m_do_zs_emulation)
      {
        // for (uint16_t sampleNum = 0; sampleNum < sam; sampleNum++)
        //   {

        for (std::unique_ptr<TpcRawHit::AdcIterator> adc_iterator(tpchit->CreateAdcIterator());
             !adc_iterator->IsDone();
             adc_iterator->Next())
        {
          // const uint16_t sampleNum = adc_iterator->CurrentTimeBin();
          const uint16_t adc = adc_iterator->CurrentAdc();

          if (adc > 0)
          {
            pedhist.Fill(adc);
          }
        }
        int hmax = 0;
        int hmaxbin = 0;
        for (int nbin = 1; nbin <= pedhist.GetNbinsX(); nbin++)
        {
          float val = pedhist.GetBinContent(nbin);
          if (val > hmax)
          {
            hmaxbin = nbin;
            hmax = val;
          }
        }

        // calculate pedestal mean and sigma

        if (pedhist.GetStdDev() == 0 || pedhist.GetEntries() == 0)
        {
          hpedestal = pedhist.GetBinCenter(pedhist.GetMaximumBin());
          hpedwidth = 999;
        }
        else
        {
          // calc peak position
          double adc_sum = 0.0;
          double ibin_sum = 0.0;
          double ibin2_sum = 0.0;

          for (int isum = -3; isum <= 3; isum++)
          {
            float val = pedhist.GetBinContent(hmaxbin + isum);
            float center = pedhist.GetBinCenter(hmaxbin + isum);
            ibin_sum += center * val;
            ibin2_sum += center * center * val;
            adc_sum += val;
          }

          hpedestal = ibin_sum / adc_sum;
          hpedwidth = sqrt(ibin2_sum / adc_sum - (hpedestal * hpedestal));
        }
        if (m_do_baseline_corr)
        {
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
          // find or insert TH2I;
          std::map<unsigned int, TH2I*>::iterator fee_map_it;

          fee_map_it = feeadc_map.find(fee_key);
          if (fee_map_it != feeadc_map.end())
          {
            feehist = (*fee_map_it).second;
          }
          else
          {
            std::string histname = "h" + std::to_string(fee_key);
            feehist = new TH2I(histname.c_str(), "histname", max_time_range + 1, -0.5, max_time_range + 0.5, 501, -0.5, 1000.5);
            feeadc_map.insert(std::make_pair(fee_key, feehist));
          }
        }
        ntotalchannels++;
        if (m_do_noise_rejection && !m_do_baseline_corr)
        {
          if (hpedwidth < 0.5 || hpedestal < 10 || hpedwidth == 999)
          {
            n_noisychannels++;
            continue;
          }
        }
      }
      else
      {
        hpedestal = 60;
        hpedwidth = m_zs_threshold;
      }
      // for (uint16_t s = 0; s < sam; s++)
      // {

      for (std::unique_ptr<TpcRawHit::AdcIterator> adc_iterator(tpchit->CreateAdcIterator());
           !adc_iterator->IsDone();
           adc_iterator->Next())
      {
        const uint16_t s = adc_iterator->CurrentTimeBin();
        const uint16_t adc = adc_iterator->CurrentAdc();
        int t = s - m_presampleShift;
        if (t < 0)
        {
          continue;
        }
        if (m_do_baseline_corr && feehist != nullptr && (!m_do_zs_emulation))
        {
          if (adc > 0)
          {
            feehist->Fill(t, adc - hpedestal + pedestal_offset);
          }
        }
        float threshold_cut = (hpedwidth * m_ped_sig_cut);
        if (m_do_zs_emulation)
        {
          threshold_cut = m_zs_threshold;
        }
        if ((float(adc) - hpedestal) > threshold_cut)
        {
          hit_key = TpcDefs::genHitKey(phibin, (unsigned int) t);
          // find existing hit, or create new one
          hit = hit_set_container_itr->second->getHit(hit_key);
          if (!hit)
          {
            hit = new TrkrHitv2();
            if (m_do_baseline_corr)
            {
              hit->setAdc(float(adc) - hpedestal + pedestal_offset);
            }
            else
            {
              hit->setAdc(float(adc) - hpedestal);
            }
            hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
          }
          if (m_writeTree)
          {
            float fXh[18];
            int nh = 0;

            fXh[nh++] = _ievent - 1;
            fXh[nh++] = 0;                        // gtm_bco;
            fXh[nh++] = 0;                        // packet_id;
            fXh[nh++] = 0;                        // ep;
            fXh[nh++] = mc_sectors[sector % 12];  // Sector;
            fXh[nh++] = side;
            fXh[nh++] = fee;
            fXh[nh++] = 0;  // channel;
            fXh[nh++] = 0;  // sampadd;
            fXh[nh++] = 0;  // sampch;
            fXh[nh++] = (float) phibin;
            fXh[nh++] = (float) t;
            fXh[nh++] = layer;
            fXh[nh++] = (float(adc) - hpedestal + pedestal_offset);
            fXh[nh++] = hpedestal;
            fXh[nh++] = hpedwidth;

            m_ntup_hits->Fill(fXh);
          }
        }
      }
    }
  }

  if (m_do_noise_rejection && Verbosity() >= 2)
  {
    std::cout << " noisy / total channels = " << n_noisychannels << "/" << ntotalchannels << " = " << n_noisychannels / (double) ntotalchannels << std::endl;
  }
  if (m_do_baseline_corr == true)
  {
    // Histos filled now process them for fee local baselines

    for (auto& hiter : feeadc_map)
    {
      if (hiter.second != nullptr)
      {
        TH2I* hist2d = hiter.second;
        std::vector<float> pedvec(hist2d->GetNbinsX(), 0);
        feebaseline_map.insert(std::make_pair(hiter.first, pedvec));
        std::map<unsigned int, std::vector<float>>::iterator fee_blm_it = feebaseline_map.find(hiter.first);
        (*fee_blm_it).second.resize(hist2d->GetNbinsX(), 0);

        for (int binx = 1; binx < hist2d->GetNbinsX(); binx++)
        {
          double timebin = ((TAxis*) hist2d->GetXaxis())->GetBinCenter(binx);
          std::string histname1d = "h" + std::to_string(hiter.first) + "_" + std::to_string((int) timebin);
          TH1D* hist1d = hist2d->ProjectionY(histname1d.c_str(), binx, binx);
          float local_ped = 0;
#ifdef DEBUG
          //  if((*hiter).first == 210802&&timebin==383){

          std::cout << " fedkey: " << (*hiter).first
                    << " entries: " << hist1d->GetEntries()
                    << std::endl;
          // }/**/
#endif

          if (hist1d->GetEntries() > 0)
          {
            int maxbin = hist1d->GetMaximumBin();
            // calc peak position
            double hadc_sum = 0.0;
            double hibin_sum = 0.0;
            //	    double hibin2_sum = 0.0;

            for (int isum = -3; isum <= 3; isum++)
            {
              float val = hist1d->GetBinContent(maxbin + isum);
              float center = hist1d->GetBinCenter(maxbin + isum);
              hibin_sum += center * val;
              // hibin2_sum += center * center * val;
              hadc_sum += val;
#ifdef DEBUG
              if ((*hiter).first == 210802 && timebin == 383)
              {
                std::cout << " fedkey: " << (*hiter).first
                          << " tbin: " << timebin
                          << " maxb " << maxbin
                          << " val: " << val
                          << " center: " << center
                          << std::endl;
              } /**/
#endif
            }
            local_ped = hibin_sum / hadc_sum;
          }
#ifdef DEBUG
          /**/
          if ((*hiter).first == 210802 && timebin == 383)
          {
            std::cout << " fedkey: " << (*hiter).first
                      << " root bin: " << binx
                      << " tbin: " << timebin
                      << " loc_ped: " << local_ped
                      << " entries: " << hist1d->GetEntries()
                      << std::endl;
          } /**/
#endif
          delete hist1d;
          (*fee_blm_it).second[(int) timebin] = local_ped;
        }
        // feebaseline_map.insert(std::make_pair((*hiter).first,pedvec));
      }
    }
    if (Verbosity() >= 1)
    {
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
        float hpedestal2 = 0;
        float hpedwidth2 = 0;
        std::map<unsigned int, chan_info>::iterator chan_it = chan_map.find(pad_key);
        if (chan_it != chan_map.end())
        {
          chan_info cinfo = (*chan_it).second;
          fee = cinfo.fee;
          hpedestal2 = cinfo.ped;
          hpedwidth2 = cinfo.width;
        }

        int rx = get_rx(layer);
        float corr = 0;

        unsigned int fee_key = create_fee_key(side, sector, rx, fee);
        std::map<unsigned int, std::vector<float>>::iterator fee_blm_it = feebaseline_map.find(fee_key);
        if (fee_blm_it != feebaseline_map.end())
        {
          corr = (*fee_blm_it).second[tbin] - pedestal_offset;
        }
        else
        {
          continue;
#ifdef DEBUG
          std::cout << " shit " << _ievent - 1
                    << " fedkey: " << fee_key
                    << " padkey: " << pad_key
                    << " layer: " << layer
                    << " side " << side
                    << " sector " << sector
                    << " fee " << fee
                    << " tbin: " << tbin
                    << " phibin " << phibin
                    << " adc " << adc
                    << std::endl;
#endif
        }
#ifdef DEBUG
        if (tbin == 383 && layer >= 7 + 32 && fee == 21)
        {
          std::cout << " before shit " << _ievent - 1
                    << " fedkey: " << fee_key
                    << " padkey: " << pad_key
                    << " layer: " << layer
                    << " side " << side
                    << " sector " << sector
                    << " fee " << fee
                    << " tbin: " << tbin
                    << " phibin " << phibin
                    << " adc " << adc
                    << std::endl;
        }
#endif
        hitr->second->setAdc(0);
        if (m_do_noise_rejection)
        {
          if (hpedwidth2 < 0.5 || hpedestal2 < 10 || hpedwidth2 == 999)
          {
            n_noisychannels++;
            continue;
          }
        }
        if (hpedwidth2 > -100 && hpedestal2 > -100)
        {
          if ((float(adc) - pedestal_offset - corr) > (hpedwidth2 * m_ped_sig_cut))
          {
            float nuadc = (float(adc) - corr - pedestal_offset);
            if (nuadc < 0)
            {
              nuadc = 0;
            }
            hitr->second->setAdc(float(nuadc));
#ifdef DEBUG
            //	    hitr->second->setAdc(10);
            if (tbin == 383 && layer >= 7 + 32 && fee == 21)
            {
              std::cout << " after shit " << _ievent - 1
                        << " fedkey: " << fee_key
                        << " padkey: " << pad_key
                        << " layer: " << layer
                        << " side " << side
                        << " sector " << sector
                        << " fee " << fee
                        << " tbin: " << tbin
                        << " phibin " << phibin
                        << " adc " << adc
                        << " corr: " << corr
                        << " adcnu " << (float(adc) - corr - pedestal_offset)
                        << " adc in " << hitr->second->getAdc()
                        << std::endl;
            }
#endif
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
              fXh[nh++] = hpedestal2;
              fXh[nh++] = hpedwidth2;
              fXh[nh++] = corr;

              m_ntup_hits_corr->Fill(fXh);
            }
          }
        }
      }
    }
  }
  // reset histogramms
  for (auto& hiter : feeadc_map)
  {
    if (hiter.second != nullptr)
    {
      hiter.second->Reset();
    }
  }
  feebaseline_map.clear();

  if (Verbosity())
  {
    std::cout << " event BCO: " << bco_min << " - " << bco_max << std::endl;
    std::cout << "TpcCombinedRawDataUnpackerDebug:: done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCombinedRawDataUnpackerDebug::End(PHCompositeNode* /*topNode*/)
{
  if (m_writeTree)
  {
    m_file->cd();
    m_ntup->Write();
    m_ntup_hits->Write();
    m_ntup_hits_corr->Write();
    m_file->Close();
  }
  if (Verbosity())
  {
    std::cout << "TpcCombinedRawDataUnpackerDebug::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  // if(m_Debug==1) hm->dumpHistos(m_filename, "RECREATE");

  return Fun4AllReturnCodes::EVENT_OK;
}
