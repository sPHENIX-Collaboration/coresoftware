#include "Tools.h"

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <mvtx/CylinderGeom_Mvtx.h>
#include <intt/CylinderGeomIntt.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <iostream>
#include <phool/phool.h>  // for PHWHERE
#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/PHObject.h>  // for PHObject
#include <trackbase_historic/SvtxTrack.h>

using std::cout;
using std::endl;

namespace G4Eval {

  //Implementation of Cluster comparator
  TrkrClusterComparer::TrkrClusterComparer (float _nphi_widths, float _nz_widths ) 
    : m_nphi_widths { _nphi_widths },
    m_nz_widths   { _nz_widths }
  {};
  int TrkrClusterComparer::init(PHCompositeNode* topNode,
      std::string name_phg4_clusters,
      std::string name_reco_clusters)
  {
    // fill bin/pixel sizes
    // ------ MVTX data ------
    PHG4CylinderGeomContainer* geom_container_mvtx 
      = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
    if (!geom_container_mvtx) {
      std::cout << PHWHERE << " Could not locate CYLINDERGEOM_MVTX " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    for (int layer=0; layer<3; ++layer) {
      auto layergeom = dynamic_cast<CylinderGeom_Mvtx *>
        (geom_container_mvtx->GetLayerGeom(layer));
      const double pitch = layergeom->get_pixel_x();
      const double length = layergeom->get_pixel_z();
      m_phistep    [layer] = pitch;
      if (layer == 0) m_zstep_mvtx         = length;
    }

    // ------ INTT data ------
    PHG4CylinderGeomContainer* geom_container_intt 
      = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
    if (!geom_container_intt) {
      std::cout << PHWHERE << " Could not locate CYLINDERGEOM_INTT " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    // get phi and Z steps for intt
    for (int layer=3; layer<7; ++layer) {
      CylinderGeomIntt* geom = 
        dynamic_cast<CylinderGeomIntt*>(geom_container_intt->GetLayerGeom(layer));
      float pitch = geom->get_strip_y_spacing();
      m_phistep [layer] = pitch;
    }

    // ------ TPC data ------
    auto geom_tpc =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    if (!geom_tpc)
    {
      std::cout << PHWHERE << " Could not locate CYLINDERCELLGEOM_SVTX node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    for (int layer=7; layer<55; ++layer) {
      PHG4TpcCylinderGeom *layergeom = geom_tpc->GetLayerCellGeom(layer);
      if (layer==7) m_zstep_tpc = layergeom->get_zstep();
      m_phistep[layer] = layergeom->get_phistep();
    }

    m_TruthClusters = 
      findNode::getClass<TrkrClusterContainer>(topNode, name_phg4_clusters.c_str());
    if (!m_TruthClusters)
    {
      std::cout << PHWHERE << " Could not locate " << name_phg4_clusters << " node" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_RecoClusters = 
      findNode::getClass<TrkrClusterContainer>(topNode, name_reco_clusters.c_str());
    if (!m_TruthClusters)
    {
      std::cout << PHWHERE << " Could not locate " << name_reco_clusters << " node" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }

  std::pair<bool, float> TrkrClusterComparer::operator()
      (TrkrDefs::cluskey key_T, TrkrDefs::cluskey key_R) 
  {
    layer = TrkrDefs::getLayer(key_T);
    if (layer > 55) {
      std::cout << " Error! Trying to compar cluster in layer > 55, "
        << "which is not programmed yet!" << std::endl;
      return {false, 0.};
    }
    
    in_mvtx = ( layer <3 );
    in_intt = (layer >2  && layer <7);
    in_tpc  = (layer >6  && layer <55);

    float phi_step = m_phistep[layer];
    float z_step = in_mvtx ? m_zstep_mvtx : m_zstep_tpc;

    clus_T = m_TruthClusters ->findCluster(key_T);
    clus_R = m_RecoClusters  ->findCluster(key_R);

    phi_T     = clus_T->getPosition(0);
    phi_R     = clus_R->getPosition(0);
    phisize_R = clus_R->getPhiSize() * phi_step;
    phisize_T = clus_T->getPhiSize() * phi_step; // only for user to get, if they want

    z_T      = clus_T->getPosition(1);
    z_R      = clus_R->getPosition(1);

    if (!in_intt) {
      zsize_R  = clus_R->getZSize() * z_step;
      zsize_T  = clus_R->getZSize() * z_step;
    }

    float dphi = fabs(phi_T-phi_R);
    while (dphi > M_PI) dphi = fabs(dphi-2*M_PI);
    z_delta   = fabs (z_T-z_R);

    float phi_stat = (m_nphi_widths * phisize_R );

    bool is_match = (phi_delta < phi_stat);
    float fit_statistic = (phi_delta / phi_stat);
    if (!in_intt) {
      float z_stat   = (m_nz_widths   * zsize_R   );
      is_match = is_match && (z_delta < z_stat);
      fit_statistic += z_delta / z_stat;
    }
    return { is_match, fit_statistic };
  }

  // Implementation of the iterable struct to get cluster keys from
  // a SvtxTrack. It is used like:
  // for (auto& cluskey : ClusKeyIter(svtx_track)) {
  //    ... // do things with cluster keys
  // }
  ClusKeyIter::ClusKeyIter(SvtxTrack* _track) :
    track {_track}
  , in_silicon { _track->get_silicon_seed()!=nullptr }
  , has_tpc    { _track->get_tpc_seed()!=nullptr }
  , no_data    { !in_silicon && !has_tpc }
  {
  }

  ClusKeyIter ClusKeyIter::begin() {
    ClusKeyIter iter0 { track }; 
    if (iter0.no_data) return iter0;
    if (iter0.in_silicon) {
      iter0.iter = track->get_silicon_seed()->begin_cluster_keys();
      iter0.iter_end_silicon = track->get_silicon_seed()->end_cluster_keys();
    } else if (has_tpc) {
      iter0.iter = track->get_tpc_seed()->begin_cluster_keys();
    } 
    return iter0;
  }

  ClusKeyIter ClusKeyIter::end() {
    ClusKeyIter iter0 { track }; 
    if (iter0.no_data) return iter0;
    if (has_tpc) {
      iter0.iter = track->get_tpc_seed()->end_cluster_keys();
    } else if (in_silicon) {
      iter0.iter = track->get_silicon_seed()->end_cluster_keys();
    } 
    return iter0;
  }

  void ClusKeyIter::operator++() {
    if (no_data) return;
    ++iter;
    if (in_silicon && has_tpc && iter == iter_end_silicon) {
      in_silicon = false;
      iter = track->get_tpc_seed()->begin_cluster_keys();
    }
  }

  bool ClusKeyIter::operator!=(const ClusKeyIter& rhs) {
    if (no_data) return false;
    return iter != rhs.iter;
  }

  TrkrDefs::cluskey ClusKeyIter::operator*() {
    return *iter;
  }
}
