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

using std::cout, std::endl;

double TpcClusterBuilder::square(double v) { return v*v; }
double TpcClusterBuilder::square(float  v) { return v*v; }

TpcClusterBuilder::TpcClusterBuilder
    ( TrkrClusterContainer*         _truth_cluster_container
    /* , TrkrTruthTrackContainer*      _truth_track_container */
    , ActsGeometry*                 _ActsGeometry
    , PHG4TpcCylinderGeomContainer* _geom_container
   ): m_clusterlist         { _truth_cluster_container }
    , truth_track_container { _truth_track_container   }
    , m_tGeometry           { _ActsGeometry            }
    , geom_container        { _geom_container          }
{}

void TpcClusterBuilder::cluster_and_reset() {

  if (!is_embedded_track) {
    m_hits->Reset();
    return;
  }

  /* int num_hitsets = 0; */
  TrkrHitSetContainer::ConstRange hitsetrange;
  hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  /* num_hitsets = std::distance(hitsetrange.first,hitsetrange.second); */

  //-----------------------------------------------------
  // from TpcClusterizer.cc ~line: 837
  //-----------------------------------------------------
  int count = 0;
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
      hitsetitr != hitsetrange.second;
      ++hitsetitr)
  {
    //if (count>0)continue;
    TrkrHitSet *hitset = hitsetitr->second;
    unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
    int side = TpcDefs::getSide(hitsetitr->first);
    unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
    PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);

    // instanciate new thread pair, at the end of thread vector
    thread_data data{};
    /* thread_pair_t& thread_pair = threads.emplace_back(); */
    data.layergeom          = layergeom;
    data.hitset             = hitset;
    /* data.rawhitset          = nullptr; */
    data.layer              = layer;
    data.pedestal           = pedestal;
    data.sector             = sector;
    data.side               = side;
    /* data.do_assoc           = do_hit_assoc; */
    /* data.do_wedge_emulation = do_wedge_emulation; */
    /* data.do_singles         = do_singles; */
    data.tGeometry          = m_tGeometry;
    data.maxHalfSizeT       = MaxClusterHalfSizeT;
    data.maxHalfSizePhi     = MaxClusterHalfSizePhi;
    data.sampa_tbias        = m_sampa_tbias;
    data.cluster_version    = cluster_version;
    data.verbosity          = 0;

    unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
    unsigned short NPhiBinsSector = NPhiBins/12;
    unsigned short NTBins = (unsigned short)layergeom->get_zbins();
    unsigned short NTBinsSide = NTBins;
    unsigned short NTBinsMin = 0;
    unsigned short PhiOffset = NPhiBinsSector * sector;
    unsigned short TOffset = NTBinsMin;

    m_tdriftmax = AdcClockPeriod * NTBins / 2.0;  
    data.m_tdriftmax = m_tdriftmax;

    data.phibins   = NPhiBinsSector;
    data.phioffset = PhiOffset;
    data.tbins     = NTBinsSide;
    data.toffset   = TOffset ;

    ProcessSectorData(&data);

    const auto hitsetkey = TpcDefs::genHitSetKey( data.layer, data.sector, data.side );      

    // should be only one cluster per hitsetkey... right?!?
    for( uint32_t index = 0; index < data.cluster_vector.size(); ++index )
    {
      // generate cluster key
      const auto ckey = TrkrDefs::genClusKey( hitsetkey, index );

      // get cluster
      auto cluster = data.cluster_vector[index];

      // insert in map
      m_clusterlist->addClusterSpecifyKey(ckey, cluster);
      if (current_track != nullptr) current_track->addCluster(ckey);
    }
    count++;
  }

  if (current_track != nullptr) {
    cout << " filling current track " << endl;
  }
  
  // clear the track
  m_hits->Reset();
}


void TpcClusterBuilder::ProcessSectorData(
    TpcClusterBuilder::thread_data* my_data) {
  const auto& _pedestal = my_data->pedestal;
  const auto& phibins   = my_data->phibins;
  const auto& phioffset = my_data->phioffset;
  const auto& tbins     = my_data->tbins ;
  const auto& toffset   = my_data->toffset ;
  const auto& layer     = my_data->layer ;
  //    int nhits = 0;
  // for convenience, create a 2D vector to store adc values in and initialize to zero
  std::vector<std::vector<unsigned short>> adcval(
      phibins, std::vector<unsigned short>(tbins, 0));
  std::multimap<unsigned short, ihit> all_hit_map;
  std::vector<ihit> hit_vect;

  int tbinmax = 498;
  int tbinmin = 0;
  if (my_data->do_wedge_emulation){
    if (layer>=7 && layer <22){
      int etacut = 249 - ((50+(layer-7))/105.5)*249;
      tbinmin = etacut;
      tbinmax -= etacut;
    }
    if (layer>=22 && layer <=48){
      int etacut = 249 - ((65+((40.5/26)*(layer-22)))/105.5)*249;
      tbinmin = etacut;
      tbinmax -= etacut;
    }
  }

  if ( my_data->hitset!=nullptr){
    TrkrHitSet *hitset = my_data->hitset;
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();

    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
        hitr != hitrangei.second;
        ++hitr){

      if ( TpcDefs::getPad(hitr->first) - phioffset < 0 ){
        //std::cout << "WARNING phibin out of range: " << TpcDefs::getPad(hitr->first) - phioffset << " | " << phibins << std::endl;
        continue;
      }
      if ( TpcDefs::getTBin(hitr->first) - toffset < 0 ){
        //std::cout << "WARNING tbin out of range: " << TpcDefs::getTBin(hitr->first) - toffset  << " | " << tbins <<std::endl;
      }
      unsigned short phibin  = TpcDefs::getPad(hitr->first) - phioffset;
      unsigned short tbin    = TpcDefs::getTBin(hitr->first) - toffset;
      unsigned short tbinorg = TpcDefs::getTBin(hitr->first);
      if (phibin>=phibins){
        //std::cout << "WARNING phibin out of range: " << phibin << " | " << phibins << std::endl;
        continue;
      }
      if (tbin>=tbins){
        //std::cout << "WARNING z bin out of range: " << tbin << " | " << tbins << std::endl;
        continue;
      }
      if (tbinorg>tbinmax||tbinorg<tbinmin)
        continue;
      float_t fadc = (hitr->second->getAdc()) - _pedestal; // proper int rounding +0.5
      unsigned short adc = 0;
      if (fadc>0) adc  = (unsigned short) fadc;
      if (phibin      >= phibins) continue;
      if (tbin        >= tbins)   continue; // tbin is unsigned int, <0 cannot happen

      if (adc>0){
        if (adc>5){
          ihit  thisHit;

          thisHit.iphi = phibin;
          thisHit.it   = tbin;
          thisHit.adc  = adc;
          thisHit.edge = 0;
          all_hit_map.insert(std::make_pair(adc, thisHit));
        }
        adcval[phibin][tbin] = (unsigned short) adc;
      }
    }
  }

  // std::cout << "done filling " << std::endl;
  while (all_hit_map.size()>0) {
    //std::cout << "all hit map size: " << all_hit_map.size() << std::endl;
    auto iter = all_hit_map.rbegin();
    if (iter == all_hit_map.rend()){
      break;
    }
    ihit hiHit = iter->second;
    int iphi = hiHit.iphi;
    int it = hiHit.it;
    //put all hits in the all_hit_map (sorted by adc)
    //start with highest adc hit
    // -> cluster around it and get vector of hits
    std::vector<ihit> ihit_list;
    int ntouch = 0;
    int nedge  =0;
    get_cluster(iphi, it, *my_data, adcval, ihit_list, ntouch, nedge );

    // -> calculate cluster parameters
    // -> add hits to truth association
    // remove hits from all_hit_map
    // repeat untill all_hit_map empty
    calc_cluster_parameter(ihit_list, *my_data, ntouch, nedge );
    remove_hits(ihit_list,all_hit_map, adcval);
    ihit_list.clear();
  }
}

void TpcClusterBuilder::get_cluster(int phibin, int tbin, const thread_data&
    my_data, const std::vector<std::vector<unsigned short>> &adcval,
    std::vector<ihit> &ihit_list, int &touch, int &edge)
{
  // search along phi at the peak in t

  int tup =0;
  int tdown =0;
  find_t_range(phibin, tbin, my_data, adcval, tdown, tup, touch, edge);
  //now we have the t extent of the cluster, go find the phi edges

  for(int it=tbin - tdown ; it<= tbin + tup; it++){
    int phiup = 0;
    int phidown = 0;
    find_phi_range(phibin, it, my_data, adcval, phidown, phiup, touch, edge);
    for (int iphi = phibin - phidown; iphi <= (phibin + phiup); iphi++){
      if(adcval[iphi][it]>0 && adcval[iphi][it]!=USHRT_MAX){
        ihit hit;
        hit.iphi = iphi;
        hit.it = it;
        hit.adc = adcval[iphi][it];
        if(touch>0){
          if((iphi == (phibin - phidown))||
              (iphi == (phibin + phiup))){
            hit.edge = 1;
          }
        }
        ihit_list.push_back(hit);
      }
    }
  }
  return;
}

void TpcClusterBuilder::calc_cluster_parameter(
    const std::vector<TpcClusterBuilder::ihit> &ihit_list, 
    TpcClusterBuilder::thread_data& my_data, 
    int ntouch, int nedge )
{

  // get z range from layer geometry
  /* these are used for rescaling the drift velocity */
  //const double z_min = -105.5;
  //const double z_max = 105.5;

  // loop over the hits in this cluster
  double t_sum = 0.0;
  double phi_sum = 0.0;
  double adc_sum = 0.0;
  double t2_sum = 0.0;
  double phi2_sum = 0.0;

  double radius = my_data.layergeom->get_radius();  // returns center of layer

  int phibinhi = -1;
  int phibinlo = 666666;
  int tbinhi = -1;
  int tbinlo = 666666;
  int clus_size = ihit_list.size();

  if (clus_size == 1) return;

  std::vector<TrkrDefs::hitkey> hitkeyvec;
  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
    double adc = iter->adc; 
    if (adc <= 0) continue;

    int iphi = iter->iphi + my_data.phioffset;
    int it   = iter->it   + my_data.toffset;
    if (iphi > phibinhi) phibinhi = iphi;
    if (iphi < phibinlo) phibinlo = iphi;
    if (it > tbinhi) tbinhi = it;
    if (it < tbinlo) tbinlo = it;

    // update phi sums
    double phi_center = my_data.layergeom->get_phicenter(iphi);
    phi_sum += phi_center * adc;
    phi2_sum += square(phi_center)*adc;

    // update t sums
    double t = my_data.layergeom->get_zcenter(it);
    t_sum += t*adc;
    t2_sum += square(t)*adc;

    adc_sum += adc;

    // capture the hitkeys for all adc values above a certain threshold
    TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, it);
    // if (adc>5)
    hitkeyvec.push_back(hitkey);
  }
  if (m_skip_noise && adc_sum < 10){
    hitkeyvec.clear();
    return;  // skip obvious noise "clusters"
  }  
  // This is the global position
  double clusphi = phi_sum / adc_sum;
  float  clusx   = radius  * cos(clusphi);
  float  clusy   = radius  * sin(clusphi);
  double clust   = t_sum   / adc_sum;
  // needed for surface identification
  double zdriftlength = clust * my_data.tGeometry->get_drift_velocity();
  // convert z drift length to z position in the TPC
  double clusz = my_data.m_tdriftmax * my_data.tGeometry->get_drift_velocity() - zdriftlength;

  if (my_data.side == 0) clusz = -clusz;

  const double phi_cov = phi2_sum/adc_sum - square(clusphi);
  /* const double t_cov = t2_sum/adc_sum - square(clust); */

  // Get the surface key to find the surface from the 
  TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey( my_data.layer,
      my_data.sector, my_data.side );      
  Acts::Vector3 global(clusx, clusy, clusz);
  TrkrDefs::subsurfkey subsurfkey = 0;

  Surface surface = my_data.tGeometry->get_tpc_surface_from_coords(
      tpcHitSetKey, global, subsurfkey);

  if (!surface)
  {
    /// If the surface can't be found, we can't track with it. So 
    /// just return and don't add the cluster to the container
    hitkeyvec.clear();
    return;
  }

  // Estimate the errors
  const double phi_err_square = (phibinhi == phibinlo) ?
    square(radius*my_data.layergeom->get_phistep())/12:
    square(radius)*phi_cov/(adc_sum*0.14);

  /* const double t_err_square = (tbinhi == tbinlo) ? */
  /*   square(my_data.layergeom->get_zstep())/12: */
  /*   t_cov/(adc_sum*0.14); */

  char tsize = tbinhi - tbinlo + 1;
  char phisize = phibinhi - phibinlo + 1;

  // SAMPA shaping bias correction
  clust = clust + my_data.sampa_tbias;

  /// convert to Acts units
  global *= Acts::UnitConstants::cm;

  Acts::Vector3 local =
    surface->transform(my_data.tGeometry->geometry().getGeoContext()).inverse()
    * global;
  local /= Acts::UnitConstants::cm;     

  // we need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)

  /* if (my_data.cluster_version==3) { */
  /*   // Fill in the cluster details */
  /*   //================ */
  /*   auto clus = new TrkrClusterv3; */
  /*   //auto clus = std::make_unique<TrkrClusterv3>(); */
  /*   clus->setAdc(adc_sum); */      
  /*   clus->setSubSurfKey(subsurfkey); */      
  /*   clus->setLocalX(local(0)); */
  /*   clus->setLocalY(clust); */
  /*   clus->setActsLocalError(0,0, phi_err_square); */
  /*   clus->setActsLocalError(1,0, 0); */
  /*   clus->setActsLocalError(0,1, 0); */
  /*   clus->setActsLocalError(1,1, t_err_square * pow(my_data.tGeometry->get_drift_velocity(),2)); */
  /*   my_data.cluster_vector.push_back(clus); */
  /* if (my_data.cluster_version==4) { */
    //	std::cout << "clus num" << my_data.cluster_vector.size() 
    //	<< " X " << local(0) << " Y " << clust << std::endl;
  if (sqrt(phi_err_square) > 0.01) {
    auto clus = new TrkrClusterv4; // 
    //auto clus = std::make_unique<TrkrClusterv3>();
    clus->setAdc(adc_sum);  
    clus->setOverlap(ntouch);
    clus->setEdge(nedge);
    clus->setPhiSize(phisize);
    clus->setZSize(tsize);
    clus->setSubSurfKey(subsurfkey);      
    clus->setLocalX(local(0));
    clus->setLocalY(clust);
    my_data.cluster_vector.push_back(clus);
  }
  /* } */

  //      if (my_data.do_assoc && my_data.clusterhitassoc){
  /* if (my_data.do_assoc) */
  /* { */
  /*   // get cluster index in vector. It is used to store associations, and build relevant cluster keys when filling the containers */
  /*   uint32_t index = my_data.cluster_vector.size()-1; */
  /*   for (unsigned int i = 0; i < hitkeyvec.size(); i++){ */
  /*     my_data.association_vector.emplace_back(index, hitkeyvec[i]); */
  /*   } */
  /* } */
  hitkeyvec.clear();
}

void TpcClusterBuilder::remove_hits(
    std::vector<ihit> &ihit_list,
    std::multimap<unsigned short, ihit> &all_hit_map,
    std::vector<std::vector<unsigned short>> &adcval)
{
  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
    unsigned short adc    = iter->adc; 
    unsigned short phibin = iter->iphi;
    unsigned short tbin   = iter->it;
    unsigned short edge   = iter->edge;
    remove_hit(adc,phibin,tbin,edge,all_hit_map,adcval);
  }
}

void TpcClusterBuilder::remove_hit(
    double adc, int phibin, int tbin, int edge, 
    std::multimap<unsigned short, ihit> &all_hit_map, 
    std::vector<std::vector<unsigned short>> &adcval)
{
  typedef std::multimap<unsigned short, ihit>::iterator hit_iterator;
  std::pair<hit_iterator, hit_iterator> iterpair = all_hit_map.equal_range(adc);
  hit_iterator it = iterpair.first;
  for (; it != iterpair.second; ++it) {
    if (it->second.iphi == phibin && it->second.it == tbin) { 
      all_hit_map.erase(it);
      break;
    }
  }
  if(edge)
    adcval[phibin][tbin] = USHRT_MAX;
  else
    adcval[phibin][tbin] = 0;
}

void TpcClusterBuilder::find_t_range(int phibin, int tbin, const thread_data&
    my_data, const std::vector<std::vector<unsigned short>> &adcval, int&
    tdown, int& tup, int &touch, int &edge) {

  const int FitRangeT= (int) my_data.maxHalfSizeT;
  const int NTBinsMax = (int) my_data.tbins;
  tup = 0;
  tdown = 0;
  for(int it=0; it< FitRangeT; it++){
    int ct = tbin + it;

    if(ct <= 0 || ct >= NTBinsMax){
      // tup = it;
      edge++;
      break; // truncate edge
    }

    if(adcval[phibin][ct] <= 0) {
      break;
    }
    if(adcval[phibin][ct] == USHRT_MAX) {
      touch++;
      break;
    }
    //check local minima and break at minimum.
    if(ct<NTBinsMax-4){//make sure we stay clear from the edge
      if(adcval[phibin][ct]+adcval[phibin][ct+1] < 
          adcval[phibin][ct+2]+adcval[phibin][ct+3]){//rising again
        tup = it+1;
        touch++;
        break;
      }
    }
    tup = it;
  }
  for(int it=0; it< FitRangeT; it++){
    int ct = tbin - it;
    if(ct <= 0 || ct >= NTBinsMax){
      //      tdown = it;
      edge++;
      break; // truncate edge
    }
    if(adcval[phibin][ct] <= 0) {
      break;
    }
    if(adcval[phibin][ct] == USHRT_MAX) {
      touch++;
      break;
    }
    if(ct>4){//make sure we stay clear from the edge
      if(adcval[phibin][ct]+adcval[phibin][ct-1] < 
          adcval[phibin][ct-2]+adcval[phibin][ct-3]){//rising again
        tdown = it+1;
        touch++;
        break;
      }
    }
    tdown = it;
  }
  return;
}
	
void TpcClusterBuilder::find_phi_range(int phibin, int tbin, const thread_data& my_data, const
    std::vector<std::vector<unsigned short>> &adcval, int& phidown, int&
    phiup, int &touch, int &edge)
{

  int FitRangePHI = (int) my_data.maxHalfSizePhi;
  int NPhiBinsMax = (int) my_data.phibins;
  phidown = 0;
  phiup = 0;
  for(int iphi=0; iphi< FitRangePHI; iphi++){
    int cphi = phibin + iphi;
    if(cphi < 0 || cphi >= NPhiBinsMax){
      // phiup = iphi;
      edge++;
      break; // truncate edge
    }

    //break when below minimum
    if(adcval[cphi][tbin] <= 0) {
      // phiup = iphi;
      break;
    }
    if(adcval[cphi][tbin] == USHRT_MAX) {
      touch++;
      break;
    }
    //check local minima and break at minimum.
    if(cphi<NPhiBinsMax-4){//make sure we stay clear from the edge
      if(adcval[cphi][tbin]+adcval[cphi+1][tbin] < 
          adcval[cphi+2][tbin]+adcval[cphi+3][tbin]){//rising again
        phiup = iphi+1;
        touch++;
        break;
      }
    }
    phiup = iphi;
  }

  for(int iphi=0; iphi< FitRangePHI; iphi++){
    int cphi = phibin - iphi;
    if(cphi < 0 || cphi >= NPhiBinsMax){
      // phidown = iphi;
      edge++;
      break; // truncate edge
    }

    if(adcval[cphi][tbin] <= 0) {
      //phidown = iphi;
      break;
    }
    if(adcval[cphi][tbin] == USHRT_MAX) {
      touch++;
      break;
    }
    if(cphi>4){//make sure we stay clear from the edge
      if(adcval[cphi][tbin]+adcval[cphi-1][tbin] < 
          adcval[cphi-2][tbin]+adcval[cphi-3][tbin]){//rising again
        phidown = iphi+1;
        touch++;
        break;
      }
    }
    phidown = iphi;
  }
  return;
}

void TpcClusterBuilder::set_current_track(TrkrTruthTrack* track) {
  current_track = track;
}

void TpcClusterBuilder::addhitset(
    TrkrDefs::hitsetkey hitsetkey, 
    TrkrDefs::hitkey hitkey, 
    float neffelectrons) 
{
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


//==================================
//----------------------------------
//==================================
//----------------------------------
//==================================
//----------------------------------
//==================================
//----------------------------------
//==================================
//----------------------------------
//==================================

