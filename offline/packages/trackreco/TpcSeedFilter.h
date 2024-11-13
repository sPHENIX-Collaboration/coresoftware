#ifndef TPCSEEDFILTER__H
#define TPCSEEDFILTER__H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/TrackSeed.h>

#include <iostream>
#include <map>
#include <vector>

// This module, if put in an anayslis chain, will filter through the TpcSeeds (tracks made
// of clusters in the TPC) and null out any tracks that don't pass the set cuts

/* class SvtxTrack; */
/* class SvtxTrackMap; */
class PHCompositeNode;
class TrackSeedContainer;

class TpcSeedFilter : public SubsysReco
{
  //--------------------------------------------------
  // Standard public interface
  //--------------------------------------------------
 public:
  int process_event(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;

  TpcSeedFilter(                                                                                  // Criteria to match a TrkrClusterContainer and track
                                                                                                  /*-----------------------------------------------------------------------------------
                                                                                                   * Input criteria for Truth Track (with nT clusters) to reco track (with nR clusters) :
                                                                                                   *  - nclus_min   : minimum number of clusters in tpc
                                                                                                   *  - min_pt
                                                                                                   *  - max_pt      : default -1 for no cut
                                                                                                   *  - min_eta     : default to no eta cut
                                                                                                   *  - max_eta     : default to no eta cut
                                                                                                   *  - clusters_cross_boundaries : require clusters in at least two
                                                                                                   *      different radial sectors.
                                                                                                   *
                                                                                                   * @Future cuts:
                                                                                                   *   - cut on cluster values? (small clusters, etc...)
                                                                                                   *   - cut on
                                                                                                   *--------------------------------------------------------*/
      unsigned int _nclus_min_ = 0, unsigned int _min_radial_sectors_ = 1, float _min_pt_ = -1.0  // no cut default
      ,
      float _max_pt_ = -1.0  // no cut default
      ,
      float _min_eta_ = -100.  // won't apply a cut
      ,
      float _max_eta_ = 100.  // won't apply a cut
      )
    : _nclus_min{_nclus_min_}
    , _min_radial_sectors{_min_radial_sectors_}
    , _must_span_sectors{_min_radial_sectors > 1}
    , _min_pt{_min_pt_}
    , _max_pt{_max_pt_}
    , _min_eta{_min_eta_}
    , _max_eta{_max_eta_}
    , _cut_on_min_pt{_min_pt != -1.0}
    , _cut_on_max_pt{_max_pt != -1.0}
    , _cut_on_min_eta{_min_eta != -100.}
    , _cut_on_max_eta{_max_eta != -100.} {};

 public:
  void set_min_radial_sectors(unsigned int nsec)
  {
    _min_radial_sectors = nsec;
    if (_min_radial_sectors > 1)
    {
      _must_span_sectors = true;
    }
  };

  void set_min_pt(float val)
  {
    _min_pt = val;
    _cut_on_min_pt = true;
  };
  void set_max_pt(float val)
  {
    _max_pt = val;
    _cut_on_max_pt = true;
  };

  void set_min_eta(float val)
  {
    _min_eta = val;
    _cut_on_min_eta = true;
  };
  void set_max_eta(float val)
  {
    _max_eta = val;
    _cut_on_max_eta = true;
  };

  void set_min_nclusters(unsigned int n) { _nclus_min = n; };

 private:
  unsigned int _nclus_min;
  unsigned int _min_radial_sectors;
  bool _must_span_sectors{false};
  float _min_pt;
  float _max_pt;
  float _min_eta;
  float _max_eta;
  bool _cut_on_min_pt{false};
  bool _cut_on_max_pt{false};
  bool _cut_on_min_eta{false};
  bool _cut_on_max_eta{false};
  TrackSeedContainer *_trackseeds = nullptr;
};

#endif
