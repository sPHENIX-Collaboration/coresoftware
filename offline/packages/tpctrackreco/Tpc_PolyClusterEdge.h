// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCPOLYCLUSTEREDGE_H
#define TPCTRACKRECO_TPCPOLYCLUSTEREDGE_H

#include "Tpc_PolyCluster.h"

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

// Number of sector edges (0, 1 or 2) a poly cluster touches, counted as in
// TpcClusterizer (sledge/sredge): a hit on the first or last pad of its
// sector. Charge from a track over the dead zone between sectors lands only
// on these pads, so such clusters sit at the edge pad, not on the track.
inline unsigned int tpc_poly_cluster_sector_edges(const Tpc_PolyCluster& cluster, PHG4TpcGeomContainer& geom)
{
  bool left = false;
  bool right = false;
  for (unsigned int i = 0; i < cluster.size_hits(); ++i)
  {
    const auto [hitsetkey, hitkey] = cluster.get_hit_index(i);
    PHG4TpcGeom* layergeom = geom.GetLayerCellGeom(TrkrDefs::getLayer(hitsetkey));
    if (!layergeom)
    {
      continue;
    }
    const int pads_per_sector = layergeom->get_phibins() / 12;
    if (pads_per_sector <= 0)
    {
      continue;
    }
    // pads are numbered around the whole layer in blocks of N per sector;
    // same local-pad convention as IdealPadMap::get_phi
    const int local_pad = static_cast<int>(TpcDefs::getPad(hitkey)) % pads_per_sector;
    left |= (local_pad == 0);
    right |= (local_pad == pads_per_sector - 1);
  }
  return (left ? 1U : 0U) + (right ? 1U : 0U);
}

#endif
