// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCSEEDSQA_H
#define TPCSEEDSQA_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

class SvtxTrack;
class PHCompositeNode;

class TpcSeedsQA : public SubsysReco
{
 public:
  TpcSeedsQA(const std::string &name = "TpcSeedsQA");

  ~TpcSeedsQA() override = default;

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

#endif  // TPCSEEDSQA_H
