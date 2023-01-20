#include "TpcClusterBuilder.h"
/* #include <trackbase/TrkrClusterv3.h> */
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TpcDefs.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <trackbase/TrkrTruthTrack.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <algorithm>

#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit
#include <cmath>  // for sqrt, cos, sin
#include <map>  // for _Rb_tree_cons...

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrTruthTrackContainer.h>

/* class TpcClusterBuilder; */

using std::cout, std::endl, std::string, std::ofstream;

double TpcClusterBuilder::square(double v) { return v*v; }
double TpcClusterBuilder::square(float  v) { return v*v; }

void TpcClusterBuilder::set_verbosity(int verbosity_level) { verbosity = verbosity_level; }

TpcClusterBuilder::TpcClusterBuilder
    ( TrkrClusterContainer*         _truth_cluster_container
    /* , TrkrTruthTrackContainer*      _truth_track_container */
    , ActsGeometry*                 _ActsGeometry
    , PHG4TpcCylinderGeomContainer* _geom_container
   ): m_clusterlist         { _truth_cluster_container }
    , m_tGeometry           { _ActsGeometry            }
    , geom_container        { _geom_container          }
{ }

void TpcClusterBuilder::cluster_and_reset(bool clear_hitsetkey_cnt) {
  if (verbosity) if (m_hits == nullptr) cout << " m_hits == nullptr! " << endl;

  if (!is_embedded_track) {
    reset(clear_hitsetkey_cnt);
    return;
  }

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
    //if (count>0)continue;
    TrkrDefs::hitsetkey hitsetkey  = hitsetitr->first;
    TrkrHitSet *hitset             = hitsetitr->second;
    unsigned int layer             = TrkrDefs::getLayer(hitsetitr->first);
    /* int side                       = TpcDefs::getSide(hitsetitr->first); */
    unsigned int sector            = TpcDefs::getSectorId(hitsetitr->first);
    PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);

    //get the maximum and minimum phi and time
    unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
    unsigned short NPhiBinsSector = NPhiBins/12;
    unsigned short NTBins = (unsigned short)layergeom->get_zbins();
    unsigned short NTBinsSide = NTBins;
    unsigned short NTBinsMin = 0;
    unsigned short PhiOffset = NPhiBinsSector * sector;
    unsigned short TOffset = NTBinsMin;

    double m_tdriftmax = AdcClockPeriod * NTBins / 2.0;  

    unsigned short phibins   = NPhiBinsSector;
    unsigned short phioffset = PhiOffset;
    unsigned short tbins     = NTBinsSide;
    unsigned short toffset   = TOffset ;
   
    // loop over the hits in this cluster
    double t_sum = 0.0;
    double phi_sum = 0.0;
    double adc_sum = 0.0;
    // double t2_sum = 0.0;
    // double phi2_sum = 0.0;

    double radius = layergeom->get_radius();  // returns center of layer

    int phibinhi = -1;
    int phibinlo = 666666;
    int tbinhi = -1;
    int tbinlo = 666666;

    auto ihit_list = hitset->getHits();
    for(auto iter = ihit_list.first; iter != ihit_list.second; ++iter){
      unsigned int adc = iter->second->getAdc(); 
      if (adc <= 0) continue;

      int iphi = TpcDefs::getPad(iter->first) - phioffset;
      int it   = TpcDefs::getTBin(iter->first) - toffset;

      if(iphi<0 || iphi>=phibins){
        //std::cout << "WARNING phibin out of range: " << phibin << " | " << phibins << std::endl;
        continue;
      }
      if(it<0 || it>=tbins){
        //std::cout << "WARNING z bin out of range: " << tbin << " | " << tbins << std::endl;
        continue;
      }

      if (iphi > phibinhi) phibinhi = iphi;
      if (iphi < phibinlo) phibinlo = iphi;
      if (it > tbinhi) tbinhi = it;
      if (it < tbinlo) tbinlo = it;

      // update phi sums
      double phi_center = layergeom->get_phicenter(iphi);
      phi_sum += phi_center * adc;
      // phi2_sum += square(phi_center)*adc;

      // update t sums
      double t = layergeom->get_zcenter(it);
      t_sum += t*adc;
      // t2_sum += square(t)*adc;

      adc_sum += adc;
    }

    // This is the global position
    double clusphi = phi_sum / adc_sum;
    float  clusx   = radius  * cos(clusphi);
    float  clusy   = radius  * sin(clusphi);
    double clust   = t_sum   / adc_sum;
    double zdriftlength = clust * m_tGeometry->get_drift_velocity();
    // convert z drift length to z position in the TPC
    double clusz = m_tdriftmax * m_tGeometry->get_drift_velocity() - zdriftlength;
    /* const double phi_cov = phi2_sum/adc_sum - square(clusphi); */

    char tsize = tbinhi - tbinlo + 1;
    char phisize = phibinhi - phibinlo + 1;

    // get the global vector3 to then get the surface local phi and z
    Acts::Vector3 global(clusx, clusy, clusz);
    TrkrDefs::subsurfkey subsurfkey = 0;

    Surface surface = m_tGeometry->get_tpc_surface_from_coords(
      hitsetkey, global, subsurfkey);

    if (!surface) {
      if (verbosity) std::cout << "Can't find the surface! with hitsetkey " << ((int)hitsetkey) << std::endl;
      continue;
    }

    global *= Acts::UnitConstants::cm;

    Acts::Vector3 local=surface->transform(m_tGeometry->geometry().getGeoContext()).inverse() * global;
    local /= Acts::UnitConstants::cm;     

    auto cluster = new TrkrClusterv4; // 
    cluster->setAdc(adc_sum);  
    /* cluster->setOverlap(ntouch); */
    /* cluster->setEdge(nedge); */
    cluster->setPhiSize(phisize);
    cluster->setZSize(tsize);
    cluster->setSubSurfKey(subsurfkey);      
    cluster->setLocalX(local(0));
    cluster->setLocalY(clust);

    // add the clusterkey
    if (hitsetkey_cnt.find(hitsetkey)==hitsetkey_cnt.end()) {
      hitsetkey_cnt[hitsetkey] = 0;
    } else {
      hitsetkey_cnt[hitsetkey] +=1;
    }
    TrkrDefs::cluskey cluskey = TrkrDefs::genClusKey(hitsetkey, hitsetkey_cnt[hitsetkey]);
    m_clusterlist->addClusterSpecifyKey(cluskey, cluster);

    if (false) { // debug print statement
        /* cout << hitsetkey_cnt[hitsetkey] << " " << cluster->getLocalX() << std::endl; */
    }

    if (current_track != nullptr) current_track->addCluster(cluskey);
  }

  reset(clear_hitsetkey_cnt);
}

void TpcClusterBuilder::set_current_track(TrkrTruthTrack* track) {
  current_track = track;
  ++ n_tracks;
}

void TpcClusterBuilder::addhitset(
    TrkrDefs::hitsetkey hitsetkey, 
    TrkrDefs::hitkey hitkey, 
    float neffelectrons) 
{
  // copy of code in PHG4TpcPadPlaneReadout::MapToPadPlane, with a switch
  // to ignore non embedded tracks
  if (!is_embedded_track) return;

  // Add the hitset to the current embedded track
  // Code from PHG4TpcPadPlaneReadout::MapToPadPlane (around lines {}.cc::386-401)
  TrkrHitSetContainer::Iterator hitsetit = m_hits->findOrAddHitSet(hitsetkey);
  // See if this hit already exists
  TrkrHit *hit = nullptr;
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

void TpcClusterBuilder::reset(bool clear_hitsetkey_cnt) {
  // clear the hitsets before the next track
  m_hits->Reset();

  // if done with the event, also clear the hitsetkey_cnt
  if (clear_hitsetkey_cnt) {
    hitsetkey_cnt.clear();
    n_tracks = 0;
  }
}

void TpcClusterBuilder::print(
    TrkrTruthTrackContainer* truth_tracks, int nclusprint) {
  cout << " ------------- content of TrkrTruthTrackContainer ---------- " << endl;
  auto& tracks = truth_tracks->getTruthTracks();
  cout << " Number of tracks: " << tracks.size() << endl;
  for (auto& track : tracks) {
    cout << " id( " << track->getTrackid() << ")  phi:eta:pt("<<
      track->getPhi()<<":"<<track->getPseudoRapidity()<<":"<<track->getPt() << ") nclusters(" 
      << track->getClusters().size() <<") ";
    int nclus = 0;
    for (auto cluskey : track->getClusters()) {
      cout << " " 
        << ((int) TrkrDefs::getHitSetKeyFromClusKey(cluskey)) <<":index(" <<
        ((int)  TrkrDefs::getClusIndex(cluskey)) << ")";
      ++nclus;
      if (nclusprint > 0 && nclus >= nclusprint) {
        cout << " ... "; 
        break;
      }
    }
    cout << endl;
  }
  cout << " ----- end of tracks in TrkrrTruthTrackContainer ------ " << endl;
}

void TpcClusterBuilder::print_file(
    TrkrTruthTrackContainer* truth_tracks, string ofile_name)
{
  ofstream fout;
  fout.open(ofile_name.c_str());
  fout << " ------------- content of TrkrTruthTrackContainer ---------- " << endl;
  auto& tracks = truth_tracks->getTruthTracks();
  fout << " Number of tracks: " << tracks.size() << endl;
  for (auto& track : tracks) {
    fout << " id( " << track->getTrackid() << ")  phi:eta:pt("<<
      track->getPhi()<<":"<<track->getPseudoRapidity()<<":"<<track->getPt() << ") nclusters(" 
      << track->getClusters().size() <<") ";
    int nclus = 0;
    for (auto cluskey : track->getClusters()) {
      auto C = m_clusterlist->findCluster(cluskey);
      fout << " " 
        << ((int) TrkrDefs::getHitSetKeyFromClusKey(cluskey)) <<":" <<
        ((int)  TrkrDefs::getClusIndex(cluskey)) << "->adc:X:phisize:Y:zsize("
        << C->getAdc()     <<":"
        << C->getLocalX()  <<":"
        << C->getPhiSize() <<":"
        << C->getLocalY()  <<":"
        << C->getZSize()  <<") ";
      ++nclus;
    }
    fout << endl;
  }
  fout << " ----- end of tracks in TrkrrTruthTrackContainer ------ " << endl;
  fout.close();
}

