// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  PHGhostRejection
 *  \brief		Class for deciding which track based on a given TPC seed is the best one
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef PHGHOSTREJECTION_H
#define PHGHOSTREJECTION_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase_historic/TrackSeed_v2.h>


#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class TrackSeedContainer;
class TrkrCluster;
class TrkrClusterContainer;

class PHGhostRejection
{
 public:
  /* PHGhostRejection() {} */
  PHGhostRejection(unsigned int verbosity, const std::vector<TrackSeed_v2>& _seeds)
    : m_verbosity { verbosity }
    , seeds { _seeds }
    , m_rejected { std::vector<bool> (seeds.size(), false) }
  {};

  // cerbosity
  void verbosity(int verb) { m_verbosity = verb; }

  // cut because too few clusters or not spanning sector boundary
  // interally updates m_rejected
  bool cut_from_clusters(int itrack);

  // cut on the ghosts: note that this also ignores all seeds failing
  // ``cut_from_clusters'' and uses the pt_cut
  void find_ghosts(const std::vector<float>& trackChi2);
  bool is_rejected(int itrack) const { return m_rejected[itrack]; };

  bool checkClusterSharing(const TrackSeed& tr1, const TrackSeed& tr2) const;

  void set_min_pt_cut(float _ptmin) { _min_pt= _ptmin; }
  void set_must_span_sectors(bool _setting) { _must_span_sectors = _setting; }
  void set_min_clusters (int _val) { _min_clusters = _val; }
  void set_phi_cut(double d) { _phi_cut = d; }
  void set_eta_cut(double d) { _eta_cut = d; }
  void set_x_cut(double d) { _x_cut = d; }
  void set_y_cut(double d) { _y_cut = d; }
  void set_z_cut(double d) { _z_cut = d; }

 private:
  unsigned int m_verbosity;
  const std::vector<TrackSeed_v2>& seeds;
  std::vector<bool> m_rejected {}; // id
  double _phi_cut = std::numeric_limits<double>::max();
  double _eta_cut = std::numeric_limits<double>::max();
  double _x_cut = std::numeric_limits<double>::max();
  double _y_cut = std::numeric_limits<double>::max();
  double _z_cut = std::numeric_limits<double>::max();

  // cuts on minimally interesting tracks --
  // here to remove noise
  double _min_pt = 0.0;
  bool   _must_span_sectors = false;
  size_t _min_clusters = 3;


  /* TrackSeedContainer *m_trackMap = nullptr; */

  /* std::map<TrkrDefs::cluskey, Acts::Vector3> m_positions; */
};

#endif  // PHGHOSTREJECTION_H
