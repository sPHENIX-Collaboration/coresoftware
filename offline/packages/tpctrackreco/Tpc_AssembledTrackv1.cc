#include "Tpc_AssembledTrackv1.h"

ClassImp(Tpc_AssembledTrackv1)

Tpc_AssembledTrackv1::Tpc_AssembledTrackv1()
{
  Reset();
}

void Tpc_AssembledTrackv1::identify(std::ostream& os) const
{
  os << "Tpc_AssembledTrackv1:"
     << " event=" << m_event
     << " track_id=" << m_track_id
     << " side=" << m_side
     << " nsegments=" << m_nsegments
     << " nrawhits=" << m_nrawhits
     << " first_layer=" << m_first_layer
     << " last_layer=" << m_last_layer
     << " first_region=" << m_first_region
     << " last_region=" << m_last_region
     << " first_sector=" << m_first_sector
     << " last_sector=" << m_last_sector
     << " phi_slope=" << m_phi_slope
     << " phi_intercept=" << m_phi_intercept
     << " tbin_slope=" << m_tbin_slope
     << " tbin_intercept=" << m_tbin_intercept
     << " chi2_phi=" << m_chi2_phi
     << " chi2_tbin=" << m_chi2_tbin
     << " ndof_phi=" << m_ndof_phi
     << " ndof_tbin=" << m_ndof_tbin
     << " vertex_valid=" << m_vertex_valid
     << " vertex_x=" << m_vertex_x
     << " vertex_y=" << m_vertex_y
     << " vertex_r=" << m_vertex_r
     << " vertex_phi=" << m_vertex_phi
     << " vertex_tbin=" << m_vertex_tbin
     << " vertex_npairs=" << m_vertex_npairs
     << " vertex_quality=" << m_vertex_quality
     << " seed_valid=" << m_seed_valid
     << " seed_pos=(" << m_seed_x << "," << m_seed_y << "," << m_seed_z << ")"
     << " seed_mom=(" << m_seed_px << "," << m_seed_py << "," << m_seed_pz << ")"
     << " source_tracks=" << m_source_track_ids.size()
     << " hit_indices=" << m_hit_indices.size()
     << std::endl;
}

void Tpc_AssembledTrackv1::Reset()
{
  m_event = 0;
  m_track_id = 0;
  m_side = 0;

  m_nsegments = 0;
  m_nblobs = 0;
  m_nrawhits = 0;
  m_first_layer = 0;
  m_last_layer = 0;
  m_first_sector = 0;
  m_last_sector = 0;
  m_first_region = 0;
  m_last_region = 0;

  m_phi_slope = 0.0;
  m_phi_intercept = 0.0;
  m_tbin_slope = 0.0;
  m_tbin_intercept = 0.0;
  m_chi2_phi = 0.0;
  m_chi2_tbin = 0.0;
  m_ndof_phi = 0;
  m_ndof_tbin = 0;

  m_vertex_valid = 0;
  m_vertex_x = 0.0;
  m_vertex_y = 0.0;
  m_vertex_r = 0.0;
  m_vertex_phi = 0.0;
  m_vertex_tbin = 0.0;
  m_vertex_npairs = 0;
  m_vertex_quality = 0.0;

  m_seed_valid = 0;
  m_seed_x = 0.0;
  m_seed_y = 0.0;
  m_seed_z = 0.0;
  m_seed_px = 0.0;
  m_seed_py = 0.0;
  m_seed_pz = 0.0;
  m_seed_cov.assign(36, 0.0);

  m_source_track_ids.clear();
  m_source_regions.clear();
  m_source_sectors.clear();
  m_hit_indices.clear();
}

int Tpc_AssembledTrackv1::isValid() const
{
  return (m_nsegments > 0 && !m_hit_indices.empty()) ? 1 : 0;
}
