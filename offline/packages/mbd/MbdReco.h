#ifndef __MBDRECO_H__
#define __MBDRECO_H__

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

class PHCompositeNode;
class MbdEvent;
class MbdPmtContainer;
class MbdVertexMap;
class MbdOut;
class MbdGeom;
class Event;
class TF1;
class TH1;

class MbdReco : public SubsysReco
{
 public:
  MbdReco(const std::string &name = "MbdReco");

  ~MbdReco() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  int createNodes(PHCompositeNode *topNode);
  int getNodes(PHCompositeNode *topNode);

  float m_tres = 0.05;
  std::unique_ptr<TF1> m_gaussian = nullptr;

  std::unique_ptr<MbdEvent>  m_mbdevent {nullptr};
  Event                     *m_event {nullptr};
  MbdOut                    *m_mbdout {nullptr};
  MbdPmtContainer           *m_mbdpmts {nullptr};
  MbdGeom                   *m_mbdgeom {nullptr};
  MbdVertexMap              *m_mbdvtxmap {nullptr};
};

#endif  // __MBDRECO_H__
