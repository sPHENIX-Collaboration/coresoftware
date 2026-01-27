// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONSEEDMERGER_H
#define PHSILICONSEEDMERGER_H

#include <fun4all/SubsysReco.h>

#include <algorithm>
#include <string>
#include <vector>

class PHCompositeNode;
class TrackSeedContainer;

class PHSiliconSeedMerger : public SubsysReco
{
 public:
  PHSiliconSeedMerger(const std::string &name = "PHSiliconSeedMerger");

  virtual ~PHSiliconSeedMerger();

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  /**
 * Set the name of the track seed container to use when retrieving silicon tracks.
 *
 * @param name Name of the TrackSeedContainer node (defaults to "SiliconTrackSeedContainer").
 */
void trackMapName(const std::string &name) { m_trackMapName = name; }
  /**
 * Set the maximum number of overlapping clusters considered during seed merging.
 * @param nclusters Maximum number of clusters that may overlap (overlap threshold).
 */
void clusterOverlap(const unsigned int nclusters) { m_clusterOverlap = nclusters; }
  /**
 * @brief Allow merging searches to include the INTT detector.
 *
 * Configure the merger to include INTT clusters in subsequent processing by disabling the MVTX-only restriction.
 */
void searchIntt() { m_mvtxOnly = false; }
  /**
 * Enable merging of silicon seed tracks during event processing.
 *
 * When enabled, the module will merge overlapping silicon seed tracks where applicable.
 */
void mergeSeeds() { m_mergeSeeds = true; }

 private:
  int getNodes(PHCompositeNode *topNode);

  TrackSeedContainer *m_siliconTracks{nullptr};
  std::string m_trackMapName{"SiliconTrackSeedContainer"};
  /**
 * Minimum number of clusters that must be shared between two silicon track seeds
 * for them to be considered overlapping.
 *
 * Defaults to 1.
 */
unsigned int m_clusterOverlap{1};
  bool m_mergeSeeds{false};
  /**
 * Restrict seed processing to the MVTX detector only.
 *
 * When set to `true`, operations that iterate or merge silicon seed tracks
 * will be limited to clusters originating from the MVTX vertex detector.
 * When `false`, clusters from other silicon detectors are included.
 */
bool m_mvtxOnly{false};
};

#endif  // PHSILICONSEEDMERGER_H


