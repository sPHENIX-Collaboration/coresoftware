// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef COSMICTRACKQA_H
#define COSMICTRACKQA_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

class SvtxTrack;
class PHCompositeNode;

class SiliconSeedsQA : public SubsysReco
{
 public:
  SiliconSeedsQA(const std::string &name = "SiliconSeedsQA");

  ~SiliconSeedsQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;

 private:
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack *track);
  void createHistos();
  std::string getHistoPrefix() const;

  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_vertexMapName = "SvtxVertexMap";
};

#endif  // COSMICTRACKQA_H
