// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKINGITERATIONCOUNTER_H
#define TRACKINGITERATIONCOUNTER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class SvtxTrackMap;
class TrackSeed;
class TrkrClusterIterationMap;

class TrackingIterationCounter : public SubsysReco
{
 public:
  TrackingIterationCounter(const std::string &name = "TrackingIterationCounter");

  ~TrackingIterationCounter() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setTrackMapName(const std::string &name) { m_trackMapName = name; }
  void iteration(const short int iter) { m_iteration = iter; }

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  void addClustersToIterationMap(TrackSeed *seed);
  short int m_iteration = 1;

  std::string m_trackMapName = "SvtxTrackMap";
  SvtxTrackMap *m_trackMap = nullptr;
  TrkrClusterIterationMap *m_iterMap = nullptr;
};

#endif  // TRACKINGITERATIONCOUNTER_H
