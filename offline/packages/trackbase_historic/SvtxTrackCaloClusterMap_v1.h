#ifndef TRACKBASEHISTORIC_SVTXTRACKCALOCLUSTERMAPV1_H
#define TRACKBASEHISTORIC_SVTXTRACKCALOCLUSTERMAPV1_H

#include "SvtxTrackCaloClusterMap.h"

#include <cstddef>   // for size_t
#include <iostream>  // for cout, ostream

class PHObject;

class SvtxTrackCaloClusterMap_v1 : public SvtxTrackCaloClusterMap
{
 public:
  SvtxTrackCaloClusterMap_v1();
  SvtxTrackCaloClusterMap_v1(const SvtxTrackCaloClusterMap_v1& map);
  SvtxTrackCaloClusterMap_v1& operator=(const SvtxTrackCaloClusterMap_v1& map);
  ~SvtxTrackCaloClusterMap_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new SvtxTrackCaloClusterMap_v1(*this); }

  bool empty() const override { return _map.empty(); }
  size_t size() const override { return _map.size(); }
  void clear() override { Reset(); }

  const std::vector<RawCluster*> get(SvtxTrack* track) const override;
  std::vector<RawCluster*> get(SvtxTrack* track) override;
  std::vector<RawCluster*> insert(SvtxTrack* track, std::vector<RawCluster*> clusters) override;
  std::vector<RawCluster*> insert(SvtxTrack* track, RawCluster* clus) override;
  size_t erase(SvtxTrack* track) override;

  ConstIter begin() const override { return _map.begin(); }
  ConstIter find(SvtxTrack* track) const override { return _map.find(track); }
  ConstIter end() const override { return _map.end(); }

  Iter begin() override { return _map.begin(); }
  Iter find(SvtxTrack* track) override { return _map.find(track); }
  Iter end() override { return _map.end(); }

 private:
  Map _map;

  ClassDefOverride(SvtxTrackCaloClusterMap_v1, 1);
};

#endif
