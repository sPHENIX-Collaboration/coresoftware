#ifndef TPC_LASEREVENTREJECTER_H
#define TPC_LASEREVENTREJECTER_H

#include <fun4all/SubsysReco.h>

class LaserEventInfo;
class PHCompositeNode;

class LaserEventRejecter : public SubsysReco
{
 public:
  LaserEventRejecter(const std::string &name = "LaserEventRejecter");
  ~LaserEventRejecter() override = default;

  int process_event(PHCompositeNode *topNode) override;

 private:
  LaserEventInfo *m_laserEventInfo = nullptr;
};

#endif
