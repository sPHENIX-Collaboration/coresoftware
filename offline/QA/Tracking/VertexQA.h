// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef VERTEXQA_H
#define VERTEXQA_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class SvtxTrack;
class PHCompositeNode;

class VertexQA : public SubsysReco
{
 public:
  VertexQA(const std::string &name = "VertexQA");

  ~VertexQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;

  void beginRun(const int run) { m_beginRun = run; }
  void endRun(const int run) { m_endRun = run; }

 private:
  void createHistos();
  // xy slope, xy int, zr slope, zr int
  std::string getHistoPrefix() const;

  int m_event = 0;
  int m_vertices = 0;
  std::string m_vertexMapName = "SvtxVertexMap";
  int m_beginRun = 25900;
  int m_endRun = 26200;
  int m_runbins = m_endRun - m_beginRun;
};

#endif  // VERTEXQA_H
