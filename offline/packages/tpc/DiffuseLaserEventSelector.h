#ifndef DiffuseLaserEventSelector_H
#define DiffuseLaserEventSelector_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class DiffuseLaserEventSelector : public SubsysReco
{
 public:
  DiffuseLaserEventSelector(const std::string& name = "DiffuseLaserEventSelector");

  ~DiffuseLaserEventSelector() override = default;

  int process_event(PHCompositeNode* topNode) override;

  void RequireTPCDiffuseLaser(bool b) { m_requireTPCDiffuseLaser = b; }
  void RequireGL1Laser(bool b) { m_requireGL1Laser = b; }
  void RejectGL1Pileup(bool b) { m_rejectGL1Pileup = b; }

 private:
  bool m_requireTPCDiffuseLaser = true;
  bool m_requireGL1Laser = false;
  bool m_rejectGL1Pileup = true;
};

#endif