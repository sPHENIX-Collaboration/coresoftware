// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef INTTXYVERTEXFINDER_H
#define INTTXYVERTEXFINDER_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class PHCompositeNode;
class InttVertexMap;
class INTTXYvtx;

class InttXYVertexFinder : public SubsysReco
{
  typedef std::array<double, 3> VertexPos;

 public:
  explicit InttXYVertexFinder(const std::string &name = "InttXYVertexFinder");

  ~InttXYVertexFinder() override;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

  void SetPeriod(const int period) { m_period = period; }

  void SetBeamCenter(const double beamx, const double beamy);
  void SetSaveHisto(const bool savehist);
  void EnableDrawHisto(const bool enable);
  void EnableQA(const bool enable);

 private:
  int createNodes(PHCompositeNode *topNode);

 private:
  INTTXYvtx *m_inttxyvtx{nullptr};
  InttVertexMap *m_inttvertexmap{nullptr};
  int m_period{1000};

  VertexPos m_vertex_quad{0, 0, -9999.};
  VertexPos m_vertex_line{0, 0, -9999.};
};

#endif  // INTTXYVERTEXFINDER_H
