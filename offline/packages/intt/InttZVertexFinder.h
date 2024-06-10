// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef INTT_INTTZVERTEXFINDER_H
#define INTT_INTTZVERTEXFINDER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class InttVertexMap;

class INTTZvtx;

class InttZVertexFinder : public SubsysReco
{
 public:
  InttZVertexFinder(const std::string &name = "InttZVertexFinder");

  ~InttZVertexFinder() override;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

  void SetBeamCenter(const double beamx, const double beamy);
  void SetOutDirectory(const std::string &outDirectory);

  void EnableQA(const bool enableQA);
  void EnableEventDisplay(const bool enableEvtDisp);

 private:
  int createNodes(PHCompositeNode *topNode);

 private:
  INTTZvtx *m_inttzvtx{nullptr};
  InttVertexMap *m_inttvertexmap{nullptr};
};

#endif  // INTT_INTTZVERTEXFINDER_H
