// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKCONTAINERCOMBINER_H
#define TRACKCONTAINERCOMBINER_H

#include <fun4all/SubsysReco.h>

#include <string>

class TrackSeedContainer;
class SvtxTrackMap;
class PHCompositeNode;

class TrackContainerCombiner : public SubsysReco
{
 public:
  TrackContainerCombiner(const std::string &name = "TrackContainerCombiner");

  ~TrackContainerCombiner() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *) override;
  int End(PHCompositeNode *) override;

  void newContainerName(const std::string &name) { m_newContainerName = name; }
  void oldContainerName(const std::string &name) { m_oldContainerName = name; }

 private:
  void mergeSeeds();
  int getNodes(PHCompositeNode *topNode);
  std::string m_newContainerName = "SiliconTrackSeedContainer";
  std::string m_oldContainerName = "SiliconTrackSeedContainerIt1";
  bool m_seedContainer = true;
  TrackSeedContainer *m_newSeedContainer = nullptr;
  TrackSeedContainer *m_oldSeedContainer = nullptr;
};

#endif  // TRACKCONTAINERCOMBINER_H
