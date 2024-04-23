#include "TpcClusterBuilder.h"

#include <trackbase/ClusHitsVerbosev1.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <g4detectors/PHG4TpcCylinderGeom.h>

#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <boost/format.hpp>

#include <algorithm>
#include <cmath>  // for sqrt, cos, sin
#include <ios>
#include <iostream>
#include <limits>
#include <map>  // for _Rb_tree_cons...
#include <set>

namespace
{
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }
}  // namespace

void TpcClusterBuilder::cluster_hits(TrkrTruthTrack* track)
{
  // now build the TrkrCluster
  TrkrHitSetContainer::ConstRange hitsetrange;
  hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);

  //-----------------------------------------------------
  // from TpcClusterizer.cc ~line: 837
  //-----------------------------------------------------

  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    // if (count>0)continue;
    TrkrDefs::hitsetkey hitsetkey = hitsetitr->first;
    TrkrHitSet* hitset = hitsetitr->second;
    unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
    int side = TpcDefs::getSide(hitsetitr->first);
    unsigned int sector = TpcDefs::getSectorId(hitsetitr->first);
    PHG4TpcCylinderGeom* layergeom = geom_container->GetLayerCellGeom(layer);

    // get the maximum and minimum phi and time
    unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
    unsigned short NPhiBinsSector = NPhiBins / 12;
    unsigned short NTBins = (unsigned short) layergeom->get_zbins();
    unsigned short NTBinsSide = NTBins;
    unsigned short NTBinsMin = 0;
    unsigned short PhiOffset = NPhiBinsSector * sector;
    unsigned short TOffset = NTBinsMin;

    double m_tdriftmax = AdcClockPeriod * NTBins / 2.0;

    unsigned short phibins = NPhiBinsSector;
    unsigned short phioffset = PhiOffset;
    unsigned short tbins = NTBinsSide;
    unsigned short toffset = TOffset;

    // loop over the hits in this cluster
    double t_sum = 0.0;
    double adc_sum = 0.0;

    /* double phi_sum = 0.0; */
    double iphi_sum = 0.0;
    // double t2_sum = 0.0;
    // double phi2_sum = 0.0;

    double radius = layergeom->get_radius();  // returns center of layer

    int phibinhi = -1;
    int phibinlo = std::numeric_limits<int>::max();
    int tbinhi = -1;
    int tbinlo = std::numeric_limits<int>::max();

    auto ihit_list = hitset->getHits();

    const int iphi_max = phioffset + phibins;
    const int it_max = toffset + tbins;

    double sum_adc{0.};  // energy = 4 * adc at this point
                         //
    // accumulate energy from all hits that are not out of range
    for (auto iter = ihit_list.first; iter != ihit_list.second; ++iter)
    {
      int iphi = TpcDefs::getPad(iter->first);  //- phioffset;
      int it = TpcDefs::getTBin(iter->first);   // - toffset;
      if (iphi < phioffset || iphi > iphi_max || it < toffset || it > it_max)
      {
        continue;
      }
      sum_adc += iter->second->getAdc();
    }
    const double threshold = sum_adc * m_pixel_thresholdrat;

    // FIXME -- see why the hits are so scattered
    std::set<int> v_iphi, v_it;                                    // FIXME
    std::map<int, unsigned int> m_iphi, m_it, m_iphiCut, m_itCut;  // FIXME
    for (auto iter = ihit_list.first; iter != ihit_list.second; ++iter)
    {
      unsigned int adc = iter->second->getAdc();
      if (adc <= 0)
      {
        continue;
      }
      int iphi = TpcDefs::getPad(iter->first);  //- phioffset;
      int it = TpcDefs::getTBin(iter->first);   // - toffset;
      if (iphi < phioffset || iphi > iphi_max)
      {
        std::cout << "WARNING phibin out of range: " << iphi << " | " << phibins << std::endl;
        continue;
      }
      if (it < toffset || it > it_max)
      {
        std::cout << "WARNING z bin out of range: " << it << " | " << tbins << std::endl;
        continue;
      }
      if (adc < threshold)
      {
        if (mClusHitsVerbose)
        {
          auto pnew = m_iphiCut.try_emplace(iphi, adc);
          if (!pnew.second)
          {
            pnew.first->second += adc;
          }

          pnew = m_itCut.try_emplace(it, adc);
          if (!pnew.second)
          {
            pnew.first->second += adc;
          }
        }
        continue;
      }

      v_iphi.insert(iphi);
      v_it.insert(it);

      // fill m_iphi
      auto pnew = m_iphi.try_emplace(iphi, adc);
      if (!pnew.second)
      {
        pnew.first->second += adc;
      }
      pnew = m_it.try_emplace(it, adc);
      if (!pnew.second)
      {
        pnew.first->second += adc;
      }

      if (iphi > phibinhi)
      {
        phibinhi = iphi;
      }
      if (iphi < phibinlo)
      {
        phibinlo = iphi;
      }
      if (it > tbinhi)
      {
        tbinhi = it;
      }
      if (it < tbinlo)
      {
        tbinlo = it;
      }

      iphi_sum += iphi * adc;
      // phi2_sum += square(phi_center)*adc;

      // update t sums
      double t = layergeom->get_zcenter(it);
      t_sum += t * adc;
      // t2_sum += square(t)*adc;

      adc_sum += adc;
    }
    if (mClusHitsVerbose)
    {
      if (verbosity > 10)
      {
        for (auto& hit : m_iphi)
        {
          std::cout << " m_phi(" << hit.first << " : " << hit.second << ") ";
        }
      }
      /* std::cout << " MELON 0 m_iphi "; */
      /* for (auto& hit : m_iphi)   std::cout  << hit.first << "|" << hit.second << " "; */
      /* std::cout << std::endl; */
      /* std::cout << " MELON 1 m_it   "; */
      /* for (auto& hit : m_it) std::cout    << hit.first << "|" << hit.second << " "; */
      /* std::cout << std::endl; */
      for (auto& hit : m_iphi)
      {
        mClusHitsVerbose->addPhiHit(hit.first, hit.second);
      }
      for (auto& hit : m_it)
      {
        mClusHitsVerbose->addZHit(hit.first, hit.second);
      }
      for (auto& hit : m_iphiCut)
      {
        mClusHitsVerbose->addPhiCutHit(hit.first, hit.second);
      }
      for (auto& hit : m_itCut)
      {
        mClusHitsVerbose->addZCutHit(hit.first, hit.second);
      }
    }

    // This is the global position
    double clusiphi = iphi_sum / adc_sum;
    double clusphi = layergeom->get_phi(clusiphi);
    /* double clusphi = phi_sum / adc_sum; */
    float clusx = radius * cos(clusphi);
    float clusy = radius * sin(clusphi);
    double clust = t_sum / adc_sum;
    double zdriftlength = clust * m_tGeometry->get_drift_velocity();
    // convert z drift length to z position in the TPC
    double clusz = m_tdriftmax * m_tGeometry->get_drift_velocity() - zdriftlength;
    if (side == 0)
    {
      clusz = -clusz;
    }

    char tsize = tbinhi - tbinlo + 1;
    char phisize = phibinhi - phibinlo + 1;

    if (tsize < 0 && verbosity > 1)
    {
      std::cout << " FIXME z4 tsize: " << ((int) tsize) << " " << tbinlo << " to " << tbinhi << std::endl;
    }

    // -------------------------------------------------
    // -------------------------------------------------
    // debug here: FIXME
    //
    /* std::cout << " FIXME phisize " << ((int) phisize) << std::endl; */
    // FIXME
    if (false)
    {  // Printing for debugging
      if ((int) phisize > 10 || (int) tsize > 8)
      {
        int _size_phi = ((int) phisize);
        int _nbins_phi = v_iphi.size();
        int _delta_phi = abs(_size_phi - _nbins_phi);
        int _size_z = ((int) tsize);
        int _nbins_z = v_it.size();
        int _delta_z = abs(_size_z - _nbins_z);

        //        TString fmt;
        std::cout << " x|" << _delta_phi << "|" << _delta_z
                  << "| new node FIXME A1  layer(" << layer
                  << ") (nset:size) phi("
                  << _nbins_phi << ":" << _size_phi << ") z("
                  << _nbins_z << ":" << _size_z << ") "
                  << "trkId(" << track->getTrackid() << ") trkpT(" << track->getPt() << ")"
                  << std::endl;
        if (phisize > 10)
        {
          std::cout << "  iphi-from-(";
          int _prev = -1;
          double tempsum = 0.;
          for (auto _ : v_iphi)
          {
            if (_prev == -1)
            {
              std::cout << _ << "): ";
            }
            else
            {
              int _diff = ((int) _ - _prev - 1);
              if (_diff != 0)
              {
                std::cout << ">" << _diff << ">";
              }
            }
            /* std::cout << std::setprecision(2) << std::fixed; */
            double _rat = (float) m_iphi[_] / (float) adc_sum;
            std::cout << boost::str(boost::format("%.2f") % _rat) << " ";
            tempsum += _rat;
            _prev = _;
          }
          if (tempsum < 0.999)
          {
            std::cout << " Z3 sumphirat: " << tempsum;
          }
          std::cout << std::endl;
        }
        if (tsize > 8)
        {
          int _prev = -1;
          double tempsum = 0.;
          std::cout << "  iz-from-(";
          for (auto _ : v_it)
          {
            if (_prev == -1)
            {
              std::cout << _ << "): ";
            }
            else
            {
              int _diff = ((int) _ - _prev - 1);
              if (_diff != 0)
              {
                std::cout << ">" << _diff << ">";
              }
            }
            /* std::cout << std::setprecision(2) << std::fixed; */
            double _rat = (float) m_it[_] / (float) adc_sum;
            std::cout << boost::str(boost::format("%.2f") % _rat) << " ";
            tempsum += _rat;
            _prev = _;
          }
          if (tempsum < 0.999)
          {
            std::cout << " Z3 sumzrat: " << tempsum;
          }
        }
        std::cout << std::endl;
      }
    }  // end debug printing

    // get the global vector3 to then get the surface local phi and z
    Acts::Vector3 global(clusx, clusy, clusz);
    TrkrDefs::subsurfkey subsurfkey = 0;

    Surface surface = m_tGeometry->get_tpc_surface_from_coords(
        hitsetkey, global, subsurfkey);

    if (!surface)
    {
      if (verbosity)
      {
        std::cout << "Can't find the surface! with hitsetkey " << ((int) hitsetkey) << std::endl;
      }
      /* if (mClusHitsVerbose) mClusHitsVerbose->Reset(); */
      continue;
    }

    // SAMPA shaping bias correction
    clust = clust + m_sampa_tbias;

    global *= Acts::UnitConstants::cm;

    Acts::Vector3 local = surface->transform(m_tGeometry->geometry().getGeoContext()).inverse() * global;
    local /= Acts::UnitConstants::cm;

    auto cluster = new TrkrClusterv4;  //
    cluster->setAdc(adc_sum);
    /* cluster->setOverlap(ntouch); */
    /* cluster->setEdge(nedge); */
    cluster->setPhiSize(phisize);
    cluster->setZSize(tsize);
    cluster->setSubSurfKey(subsurfkey);
    cluster->setLocalX(local(0));
    cluster->setLocalY(clust);

    // add the clusterkey
    auto empl = hitsetkey_cnt.try_emplace(hitsetkey, 0);
    if (!empl.second)
    {
      empl.first->second += 1;
    }
    TrkrDefs::cluskey cluskey = TrkrDefs::genClusKey(hitsetkey, hitsetkey_cnt[hitsetkey]);
    m_clusterlist->addClusterSpecifyKey(cluskey, cluster);

    track->addCluster(cluskey);
    if (mClusHitsVerbose)
    {
      mClusHitsVerbose->push_hits(cluskey);
      /* std::cout << " FIXME z1 ClusHitsVerbose.size: " << mClusHitsVerbose->getMap().size() << std::endl; */
    }
  }
  m_hits->Reset();
}

/* void TpcClusterBuilder::set_current_track(TrkrTruthTrack* track) { */
/* current_track = track; */
/* ++ n_tracks; */
/* } */

void TpcClusterBuilder::addhitset(
    TrkrDefs::hitsetkey hitsetkey,
    TrkrDefs::hitkey hitkey,
    float neffelectrons)
{
  // copy of code in PHG4TpcPadPlaneReadout::MapToPadPlane, with a switch
  // to ignore non embedded tracks
  if (!b_collect_hits)
  {
    return;
  }

  // Add the hitset to the current embedded track
  // Code from PHG4TpcPadPlaneReadout::MapToPadPlane (around lines {}.cc::386-401)
  TrkrHitSetContainer::Iterator hitsetit = m_hits->findOrAddHitSet(hitsetkey);
  // See if this hit already exists
  TrkrHit* hit = nullptr;
  hit = hitsetit->second->getHit(hitkey);
  if (!hit)
  {
    // create a new one
    hit = new TrkrHitv2();
    hitsetit->second->addHitSpecificKey(hitkey, hit);
  }
  // Either way, add the energy to it  -- adc values will be added at digitization
  hit->addEnergy(neffelectrons);
}

void TpcClusterBuilder::clear_hitsetkey_cnt()
{
  hitsetkey_cnt.clear();
}

void TpcClusterBuilder::print(
    TrkrTruthTrackContainer* truth_tracks, int nclusprint)
{
  std::cout << " ------------- content of TrkrTruthTrackContainer ---------- " << std::endl;
  auto& tmap = truth_tracks->getMap();
  std::cout << " Number of tracks:  xyz db : " << tmap.size() << std::endl;
  for (auto& _pair : tmap)
  {
    auto& track = _pair.second;
    std::cout << boost::str(boost::format("id(%2i) phi:eta:pt(") % ((int) track->getTrackid()))
              << boost::str(boost::format("%5.2f:%5.2f:%5.2f") % track->getPhi() % track->getPseudoRapidity() % track->getPt())
              << ") nclusters(" << track->getClusters().size() << ") ";
    if (verbosity <= 10)
    {
      std::cout << std::endl;
    }
    else
    {
      int nclus = 0;
      for (auto cluskey : track->getClusters())
      {
        std::cout << " "
                  << ((int) TrkrDefs::getHitSetKeyFromClusKey(cluskey)) << ":index(" << ((int) TrkrDefs::getClusIndex(cluskey)) << ")";
        ++nclus;
        if (nclusprint > 0 && nclus >= nclusprint)
        {
          std::cout << " ... ";
          break;
        }
      }
    }
  }
  std::cout << " ----- end of tracks in TrkrrTruthTrackContainer ------ " << std::endl;
}

void TpcClusterBuilder::print_file(
    TrkrTruthTrackContainer* truth_tracks, const std::string& ofile_name)
{
  std::ofstream fout;
  fout.open(ofile_name.c_str());
  fout << " ------------- content of TrkrTruthTrackContainer ---------- " << std::endl;
  auto& tmap = truth_tracks->getMap();
  fout << " Number of tracks: " << tmap.size() << std::endl;
  for (auto& _pair : tmap)
  {
    auto& track = _pair.second;
    fout << " id( " << track->getTrackid() << ")  phi:eta:pt(" << track->getPhi() << ":" << track->getPseudoRapidity() << ":" << track->getPt() << ") nclusters("
         << track->getClusters().size() << ") ";
    int nclus = 0;
    for (auto cluskey : track->getClusters())
    {
      auto C = m_clusterlist->findCluster(cluskey);
      fout << " "
           << ((int) TrkrDefs::getHitSetKeyFromClusKey(cluskey)) << ":" << ((int) TrkrDefs::getClusIndex(cluskey)) << "->adc:X:phisize:Y:zsize("
           << C->getAdc() << ":"
           << C->getLocalX() << ":"
           << C->getPhiSize() << ":"
           << C->getLocalY() << ":"
           << C->getZSize() << ") ";
      ++nclus;
    }
    fout << std::endl;
  }
  fout << " ----- end of tracks in TrkrrTruthTrackContainer ------ " << std::endl;
  fout.close();
}

void TpcClusterBuilder::set_input_nodes(
    TrkrClusterContainer* _truth_cluster_container, ActsGeometry* ActsGeometry, PHG4TpcCylinderGeomContainer* _geom_container, ClusHitsVerbosev1* _clushitsverbose)
{
  m_clusterlist = _truth_cluster_container;
  m_tGeometry = ActsGeometry;
  geom_container = _geom_container;
  mClusHitsVerbose = _clushitsverbose;
  m_needs_input_nodes = false;
};
