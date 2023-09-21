#ifndef INTT_INTTVERTEXFINDER_H
#define INTT_INTTVERTEXFINDER_H

#include <fun4all/SubsysReco.h>

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
  double calculateZvertex(double* zcenter=nullptr, double* zrms=nullptr, double* zmean=nullptr);

 private:
  // node tree storage pointers
  InttVertexMap        *m_inttvertexmap = nullptr;
  ActsGeometry         *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr; 

  TH1*                  h_zvtxseed_ = nullptr;

  // settings
  //
  double xbeam_ = 0.;
  double ybeam_ = 0.;
};

#endif
