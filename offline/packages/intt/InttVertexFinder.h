#ifndef INTT_INTTVERTEXFINDER_H
#define INTT_INTTVERTEXFINDER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>

#include <string>

class PHCompositeNode;
class ActsGeometry;
class TrkrClusterContainer;
class TH1;
class InttVertexMap;

class InttVertexFinder : public SubsysReco
{
 public:
  InttVertexFinder(const std::string &name = "InttVertexFinder");
  ~InttVertexFinder() override;

  //! module initialization
  int Init(PHCompositeNode */*topNode*/) override;

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! end of process
  int End(PHCompositeNode */*topNode*/) override { return 0; }


  void setBeamCenter(double x=0, double y=0) {xbeam_=x; ybeam_=y;}


 private:
  int    createNodes(PHCompositeNode *topNode);
  double calculateZvertex(double* zcenter=NULL, double* zrms=NULL, double* zmean=NULL);

 private:
  // node tree storage pointers
  InttVertexMap        *m_inttvertexmap;
  ActsGeometry         *m_tGeometry;
  TrkrClusterContainer *m_clusterlist; 

  TH1*                  h_zvtxseed_;

  // settings
  //
  double xbeam_, ybeam_;
};

#endif
