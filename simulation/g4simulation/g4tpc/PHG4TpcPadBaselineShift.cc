#include "PHG4TpcPadBaselineShift.h"

#include <trackbase/TpcDefs.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <trackbase/ActsSurfaceMaps.h>       // for ActsSurfaceMaps
#include <trackbase/ActsTrackingGeometry.h>  // for ActsTrackingG...
#pragma GCC diagnostic pop

#include <trackbase/TrkrClusterContainer.h>  // for TrkrClusterCo...
#include <trackbase/TrkrClusterHitAssoc.h>   // for TrkrClusterHi...
#include <trackbase/TrkrHit.h>               // for TrkrHit

#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

#include <cmath>  // for sqrt, cos, sin
#include <iostream>
#include <map>  // for _Rb_tree_cons...
#include <string>
#include <utility>  // for pair
#include <vector>

namespace
{
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }
}  // namespace

int findRBin(float R)
{
  // Finding pad number from the center (bin) for hits
  int binR = -1;
  //Realistic binning
  //double r_bins[r_bins_N+1] = {217.83,
  //                            311.05,317.92,323.31,329.27,334.63,340.59,345.95,351.91,357.27,363.23,368.59,374.55,379.91,385.87,391.23,397.19,402.49,
  //                            411.53,421.70,431.90,442.11,452.32,462.52,472.73,482.94,493.14,503.35,513.56,523.76,533.97,544.18,554.39,564.59,574.76,
  //                            583.67,594.59,605.57,616.54,627.51,638.48,649.45,660.42,671.39,682.36,693.33,704.30,715.27,726.24,737.21,748.18,759.11};
  const int r_bins_N = 53;
  double r_bins[r_bins_N + 1];
  r_bins[0] = 30.3125;
  double bin_width = 0.625;
  for (int i = 1; i < r_bins_N; i++)
  {
    if (i == 16) bin_width = 0.9375;
    if (i > 16) bin_width = 1.25;
    if (i == 31) bin_width = 1.1562;
    if (i > 31) bin_width = 1.0624;

    r_bins[i] = r_bins[i - 1] + bin_width;
  }

  double R_min = 30;
  while (R > R_min)
  {
    binR += 1;
    R_min = r_bins[binR];
  }
  return binR;
}

//____________________________________________________________________________..
PHG4TpcPadBaselineShift::PHG4TpcPadBaselineShift(const std::string &name)
  : SubsysReco(name)
{
  std::cout << "PHG4TpcPadBaselineShift::PHG4TpcPadBaselineShift(const std::string &name) Calling ctor" << std::endl;
}

bool PHG4TpcPadBaselineShift::is_in_sector_boundary(int phibin, int sector, PHG4TpcCylinderGeom *layergeom)
{
  bool reject_it = false;

  // sector boundaries occur every 1/12 of the full phi bin range
  int PhiBins = layergeom->get_phibins();
  int PhiBinsSector = PhiBins / 12;

  double radius = layergeom->get_radius();
  double PhiBinSize = 2.0 * radius * M_PI / (double) PhiBins;

  // sector starts where?
  int sector_lo = sector * PhiBinsSector;
  int sector_hi = sector_lo + PhiBinsSector - 1;

  int sector_fiducial_bins = (int) (SectorFiducialCut / PhiBinSize);

  if (phibin < sector_lo + sector_fiducial_bins || phibin > sector_hi - sector_fiducial_bins)
  {
    reject_it = true;
  }

  return reject_it;
}
//____________________________________________________________________________..
PHG4TpcPadBaselineShift::~PHG4TpcPadBaselineShift()
{
  std::cout << "PHG4TpcPadBaselineShift::~PHG4TpcPadBaselineShift() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int PHG4TpcPadBaselineShift::Init(PHCompositeNode * /*topNode*/)
{
  //outfile = new TFile(_filename.c_str(), "RECREATE");
  _hit_z = 0;
  _hit_r = 0;
  _hit_phi = 0;
  _hit_e = 0;
  _hit_adc = 0;
  _hit_adc_bls = 0;
  _hit_layer = -1;
  _hit_sector = -1;
  if (_writeTree == 1)
  {
    outfile = new TFile(_filename.c_str(), "RECREATE");
    _rawHits = new TTree("hTree", "tpc hit tree for base-line shift tests");
    _rawHits->Branch("z", &_hit_z);
    _rawHits->Branch("r", &_hit_r);
    _rawHits->Branch("phi", &_hit_phi);
    _rawHits->Branch("e", &_hit_e);
    _rawHits->Branch("adc", &_hit_adc);
    _rawHits->Branch("adc_BLS", &_hit_adc_bls);
    _rawHits->Branch("hit_layer", &_hit_layer);
    _rawHits->Branch("_hit_sector", &_hit_sector);
  }
  return 0;
  //std::cout << "PHG4TpcPadBaselineShift::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  //return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHG4TpcPadBaselineShift::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::cout << "PHG4TpcPadBaselineShift::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHG4TpcPadBaselineShift::process_event(PHCompositeNode *topNode)
{
  //  int print_layer = 18;

  if (Verbosity() > 1000)
  {
    std::cout << "PHG4TpcPadBaselineShift::Process_Event" << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster hit associations
  m_clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterhitassoc)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERHITASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHG4TpcCylinderGeomContainer *geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,
                                                         "ActsTrackingGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE
              << "ActsTrackingGeometry not found on node tree. Exiting"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode,
                                                   "ActsSurfaceMaps");
  if (!m_surfMaps)
  {
    std::cout << PHWHERE
              << "ActsSurfaceMaps not found on node tree. Exiting"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // The hits are stored in hitsets, where each hitset contains all hits in a given TPC readout (layer, sector, side), so clusters are confined to a hitset
  // The TPC clustering is more complicated than for the silicon, because we have to deal with overlapping clusters

  // loop over the TPC HitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  //const int num_hitsets = std::distance(hitsetrange.first,hitsetrange.second);

  //Loop over R positions & sectors
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    TrkrHitSet *hitset = hitsetitr->second;
    unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
    int side = TpcDefs::getSide(hitsetitr->first);
    unsigned int sector = TpcDefs::getSectorId(hitsetitr->first);
    PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);

    _hit_sector = sector;
    _hit_layer = layer;

    double radius = layergeom->get_radius();  //in cm

    _hit_r = radius;

    unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
    unsigned short NPhiBinsSector = NPhiBins / 12;
    unsigned short NZBins = (unsigned short) layergeom->get_zbins();
    unsigned short NZBinsSide = NZBins / 2;
    unsigned short NZBinsMin = 0;
    unsigned short PhiOffset = NPhiBinsSector * sector;

    if (side == 0)
    {
      NZBinsMin = 0;
      NZBinsMax = NZBins / 2 - 1;
    }
    else
    {
      NZBinsMin = NZBins / 2;
      NZBinsMax = NZBins;
    }
    unsigned short ZOffset = NZBinsMin;
    //Gives per pad ADC for particular R and sector in event
    int perPadADC = 0;
    unsigned short phibins = NPhiBinsSector;
    unsigned short phioffset = PhiOffset;
    unsigned short zbins = NZBinsSide;
    unsigned short zoffset = ZOffset;
    float sumADC = perPadADC;  // 14 lines later just set to zero

    //phibins - number of pads in the sector
    //The Sampa clock time is 18.8 MHz, so the sampling time is 53.2 ns = 1 Z bin.
    //The Z bin range for one side of the TPC is 0-248.
    //Check: 249 bins x 53.2 ns = 13.2 microseconds.
    //The maximum drift time in the TPC is 13.2 microseconds.
    //So, 53.2 ns / Z bin.

    TF1 *f1 = new TF1("f1", "[0]*exp(-(x-[1])/[2])", 0, 1000);
    f1->SetParameter(0, 0.005);
    f1->SetParameter(1, 0);
    f1->SetParameter(2, 60);  // in terms of 50nsec time bins

    //    sumADC=0.; // this is set to zero 14 lines up (perPadADC is zero)
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();

    std::vector<unsigned short> adcval(zbins, 0);
    //    std::multimap<unsigned short, ihit> all_hit_map;
    //    std::vector<ihit> hit_vect;
    // Loop over phi & z
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first; hitr != hitrangei.second; ++hitr)
    {
      unsigned short phibin = TpcDefs::getPad(hitr->first) - phioffset;
      unsigned short zbin = TpcDefs::getTBin(hitr->first) - zoffset;
      float_t fadc = (hitr->second->getAdc());  // proper int rounding +0.5
      unsigned short adc = 0;
      if (fadc > 0)
      {
        adc = (unsigned short) fadc;
      }
      if (phibin >= phibins) continue;
      if (zbin >= zbins) continue;  // zbin is unsigned int, <0 cannot happen
      adcval[zbin] = (unsigned short) adc;
      sumADC += adc;
    }
    //Define ion-induced charge
    sumADC /= phibins;
    float ind_charge = -0.5 * sumADC * _CScale;  //CScale is the coefficient related to the capacitance of the bottom layer of the bottom GEM
    double pi = M_PI;

    for (TrkrHitSet::ConstIterator hitr = hitrangei.first; hitr != hitrangei.second; ++hitr)
    {
      unsigned short phibin = TpcDefs::getPad(hitr->first);
      unsigned short tbin = TpcDefs::getTBin(hitr->first);
      // Get the hitkey
      TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(phibin, tbin);
      TrkrHit *hit = nullptr;
      hit = hitsetitr->second->getHit(hitkey);

      tbin = TpcDefs::getTBin(hitr->first);
      phibin = TpcDefs::getPad(hitr->first);
      double phi_center = layergeom->get_phicenter(phibin);
      if (phi_center < 0) phi_center += 2 * pi;
      _hit_phi = phi_center;
      _hit_z = AdcClockPeriod * MaxTBins * _drift_velocity / 2.0 - layergeom->get_zcenter(tbin) * _drift_velocity;
      if(side == 0) _hit_z *= -1.0;

      if (hit) _hit_e = hit->getEnergy();
      _hit_adc = 0;
      _hit_adc_bls = 0;

      float_t fadc = (hitr->second->getAdc());  // proper int rounding +0.5
      _hit_adc = fadc;
      _hit_adc_bls = _hit_adc + int(ind_charge);

      if (hit && _hit_adc > 0)
      {
        //Trkr hit has only one value m_adc which is energy and ADC at the same time
        if (_hit_adc_bls > 0)
        {
          hit->setAdc(_hit_adc_bls);
        }
        else
        {
          hit->setAdc(0);
        }
      }
      if (_writeTree == 1) _rawHits->Fill();
    }

    //hitsetitr++;
  }

  //pthread_attr_destroy(&attr);

  // wait for completion of all threads
  //for( const auto& thread_pair:threads )
  //{
  //  int rc2 = pthread_join(thread_pair.thread, nullptr);
  //  if (rc2)
  //  { std::cout << "Error:unable to join," << rc2 << std::endl; }
  //}

  if (Verbosity() > 0)
    std::cout << "TPC Clusterizer found " << m_clusterlist->size() << " Clusters " << std::endl;
  std::cout << "PHG4TpcPadBaselineShift::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
//int PHG4TpcPadBaselineShift::ResetEvent(PHCompositeNode *topNode)
//{
//  std::cout << "PHG4TpcPadBaselineShift::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
//int PHG4TpcPadBaselineShift::EndRun(const int runnumber)
//{
//  std::cout << "PHG4TpcPadBaselineShift::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
int PHG4TpcPadBaselineShift::End(PHCompositeNode * /*topNode*/)
{
  if (_writeTree == 1)
  {
    outfile->cd();
    outfile->Write();
    outfile->Close();
    std::cout << "PHG4TpcPadBaselineShift::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  //return 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
//int PHG4TpcPadBaselineShift::Reset(PHCompositeNode *topNode)
//{
// std::cout << "PHG4TpcPadBaselineShift::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
//void PHG4TpcPadBaselineShift::Print(const std::string &what) const
//{
//  std::cout << "PHG4TpcPadBaselineShift::Print(const std::string &what) const Printing info for " << what << std::endl;
//}
void PHG4TpcPadBaselineShift::setScale(float CScale)
{
  _CScale = CScale;
  std::cout << "PHG4TpcPadBaselineShift::setFileName: Scale factor is set to:" << CScale << std::endl;
}
void PHG4TpcPadBaselineShift::setFileName(const std::string &filename)
{
  _filename = filename;
  std::cout << "PHG4TpcPadBaselineShift::setFileName: Output file name for PHG4TpcPadBaselineShift is set to:" << filename << std::endl;
}
void PHG4TpcPadBaselineShift::writeTree(int f_writeTree)
{
  _writeTree = f_writeTree;
  std::cout << "PHG4TpcPadBaselineShift::writeTree: True" << std::endl;
}
