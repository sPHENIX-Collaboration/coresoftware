#include "TrkrClusterIsMatcher.h"
#include "g4evalfn.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <intt/CylinderGeomIntt.h>
#include <mvtx/CylinderGeom_Mvtx.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

int TrkrClusterIsMatcher::init(
    PHCompositeNode* topNode,
    const std::string& name_phg4_clusters,
    const std::string& name_reco_clusters)
{
  // ------ MVTX data ------
  PHG4CylinderGeomContainer* geom_container_mvtx = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!geom_container_mvtx)
  {
    std::cout << PHWHERE << " Could not locate CYLINDERGEOM_MVTX " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  for (int i = 0; i < 3; ++i)
  {
    auto layergeom = dynamic_cast<CylinderGeom_Mvtx*>(geom_container_mvtx->GetLayerGeom(i));
    pitch_phi[i] = layergeom->get_pixel_x();
    if (i == 0)
    {
      pitch_z_MVTX = layergeom->get_pixel_z();
    }
  }

  // ------ INTT data ------
  PHG4CylinderGeomContainer* geom_container_intt = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container_intt)
  {
    std::cout << PHWHERE << " Could not locate CYLINDERGEOM_INTT " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // get phi and Z steps for intt
  for (int i = 3; i < 7; ++i)
  {
    CylinderGeomIntt* geom =
        dynamic_cast<CylinderGeomIntt*>(geom_container_intt->GetLayerGeom(i));
    pitch_phi[i] = geom->get_strip_y_spacing();
  }

  // ------ TPC data ------
  auto geom_tpc =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_tpc)
  {
    std::cout << PHWHERE << " Could not locate CYLINDERCELLGEOM_SVTX node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  for (int i = 7; i < 55; ++i)
  {
    PHG4TpcCylinderGeom* layergeom = geom_tpc->GetLayerCellGeom(i);
    if (i == 7)
    {
      step_t_TPC = layergeom->get_zstep();
    }
    pitch_phi[i] = layergeom->get_phistep() * layergeom->get_radius();
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

  m_ActsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_ActsGeometry)
  {
    std::cout << PHWHERE << " Could not locate ActsGeometry node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  set_tol_z_MVTX(tol_z_MVTX);
  set_tol_phi_MVTX(tol_phi_MVTX);

  set_tol_phi_INTT(tol_phi_INTT);

  set_tol_t_TPC(tol_t_TPC);
  set_tol_phi_TPC(tol_phi_TPC);

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrkrClusterIsMatcher::set_tol_phi_MVTX(float _val)
{
  tol_phi_MVTX = _val;
  for (int i = 0; i < 3; ++i)
  {
    tol_pitch_phi[i] = tol_phi_MVTX * pitch_phi[i];
  }
}

void TrkrClusterIsMatcher::set_tol_z_MVTX(float _val)
{
  tol_z_MVTX = _val;
  tol_pitch_z_MVTX = tol_z_MVTX * pitch_z_MVTX;
}

void TrkrClusterIsMatcher::set_tol_phi_INTT(float _val)
{
  tol_phi_INTT = _val;
  for (int i = 3; i < 7; ++i)
  {
    tol_pitch_phi[i] = tol_phi_INTT * pitch_phi[i];
  }
}

void TrkrClusterIsMatcher::set_tol_phi_TPC(float _val)
{
  tol_phi_TPC = _val;
  for (int i = 7; i < 55; ++i)
  {
    tol_pitch_phi[i] = tol_phi_TPC * pitch_phi[i];
  }
}

void TrkrClusterIsMatcher::set_tol_t_TPC(float _val)
{
  tol_t_TPC = _val;
  tol_step_t_TPC = tol_t_TPC * step_t_TPC;
}

bool TrkrClusterIsMatcher::operator()(TrkrDefs::cluskey key_T, TrkrDefs::cluskey key_R)
{
  // note: can use returned values, or just pull from these values
  layer = TrkrDefs::getLayer(key_T);
  if (layer > 55)
  {
    std::cout << " Error! Trying to compar cluster in layer > 55, "
              << "which is not programmed yet!" << std::endl;
    return false;
  }

  clus_T = m_TruthClusters->findCluster(key_T);
  clus_R = m_RecoClusters->findCluster(key_R);

  det_0123 = g4evalfn::trklayer_det(layer);
  dphi = g4evalfn::abs_dphi(clus_T->getPosition(0), clus_R->getPosition(0));
  switch (det_0123)
  {
  case g4evalfn::DET::MVTX:
  {
    if (single_pixel_phi_MVTX)
    {
      if (dphi > tol_pitch_phi[layer])
      {
        return false;
      }
    }
    else
    {
      if (dphi > tol_pitch_phi[layer] * std::max(clus_T->getPhiSize(), clus_R->getPhiSize()))
      {
        return false;
      }
    }
    const float delta_z = fabs(clus_T->getPosition(1) - clus_R->getPosition(1));
    if (single_pixel_z_MVTX)
    {
      if (delta_z > tol_pitch_z_MVTX)
      {
        return false;
      }
    }
    else
    {
      if (delta_z > tol_pitch_z_MVTX * std::max(clus_T->getZSize(), clus_R->getZSize()))
      {
        return false;
      }
    }
    return true;
  }
  break;

  case g4evalfn::DET::INTT:
  {
    if (single_pixel_phi_INTT)
    {
      if (dphi > tol_pitch_phi[layer])
      {
        return false;
      }
    }
    else
    {
      if (dphi > tol_pitch_phi[layer] * std::max(clus_T->getPhiSize(), clus_R->getPhiSize()))
      {
        return false;
      }
    }
    return true;
  }
  break;

  case g4evalfn::DET::TPC:
  {
    if (single_bin_phi_TPC)
    {
      if (dphi > tol_pitch_phi[layer])
      {
        return false;
      }
    }
    else
    {
      if (dphi > tol_pitch_phi[layer] * std::max(clus_T->getPhiSize(), clus_R->getPhiSize()))
      {
        return false;
      }
    }
    const float delta_t = fabs(clus_T->getPosition(1) - clus_R->getPosition(1));
    if (single_bin_t_TPC)
    {
      if (delta_t > tol_step_t_TPC)
      {
        return false;
      }
    }
    else
    {
      if (delta_t > tol_step_t_TPC * std::max(clus_T->getZSize(), clus_R->getZSize()))
      {
        return false;
      }
    }
    return true;
  }
  break;

  case g4evalfn::DET::TPOT:
  {                // TPOT
    return false;  // no info for matching TPOT at this time
  }
  break;
  }
  return false;  // code shouldn't arrive here; just for completeness for compiler
}
