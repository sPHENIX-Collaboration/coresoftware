#ifndef TRACKRECO_PHTPCRESIDUALS_H
#define TRACKRECO_PHTPCRESIDUALS_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include "ActsTrack.h"
#include "ActsTrackingGeometry.h"

class PHCompositeNode;

#include <memory>
#include <map>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>

class PHTpcResiduals : public SubsysReco
{

 public:

  PHTpcResiduals(const std::string &name = "PHTpcResiduals");
  ~PHTpcResiduals();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:

  int getNodes(PHCompositeNode *topNode);

  int getTpcResiduals(PHCompositeNode *topNode);

  void calculateTpcResiduals(const std::vector<SourceLink> sourceLinks,
			     const Acts::Vector3D momentum);

  std::map<unsigned int, ActsTrack> *m_actsProtoTracks = nullptr;
  ActsTrackingGeometry *m_tGeometry;

  TFile *outfile;
  TH2 *h_rphiResid;
  TH2 *h_zResid;
  TH2 *h_etaResidLayer;
  TH2 *h_zResidLayer;
  TH2 *h_etaResid;

};


#endif
