#ifndef PHG4CHIMNEYSUBSYSTEM_H
#define PHG4CHIMNEYSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>

class PHCompositeNode;
class PHG4ChimneyDetector;
class PHG4Detector;
class PHG4DisplayAction;

class PHG4ChimneySubsystem : public PHG4DetectorSubsystem
{
 public:
  PHG4ChimneySubsystem(const std::string &name = "CHIMNEY", const int lyr = 0);
  ~PHG4ChimneySubsystem() override;

  int InitRunSubsystem(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override { return 0; }
  void Print(const std::string &what = "ALL") const override;

  PHG4Detector *GetDetector() const override;
  PHG4DisplayAction *GetDisplayAction() const override { return m_DisplayAction; }

 private:
  void SetDefaultParameters() override;

  PHG4ChimneyDetector *m_Detector = nullptr;
  PHG4DisplayAction *m_DisplayAction = nullptr;
};

#endif
