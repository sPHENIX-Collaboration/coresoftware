#ifndef TPC_LASEREVENTINFOV2_H
#define TPC_LASEREVENTINFOV2_H

#include "LaserEventInfov1.h"

#include <iostream>

class LaserEventInfov2 : public LaserEventInfov1
{
 public:
  LaserEventInfov2() = default;
  ~LaserEventInfov2();

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override;

  PHObject* CloneMe() const override { return new LaserEventInfov2(*this); }
  void CopyTo(LaserEventInfo *info) override;

  bool isGl1LaserEvent() const override { return m_isGl1LaserEvent; }
  void setIsGl1LaserEvent(const bool isGl1LaserEvent) override { m_isGl1LaserEvent = isGl1LaserEvent; }

  bool isGl1LaserPileupEvent() const override { return m_isGl1LaserPileupEvent; }
  void setIsGl1LaserPileupEvent(const bool isGl1LaserPileupEvent) override { m_isGl1LaserPileupEvent = isGl1LaserPileupEvent; }

 protected:

  bool m_isGl1LaserEvent{false};
  bool m_isGl1LaserPileupEvent{false};

  ClassDefOverride(LaserEventInfov2, 1);
};

#endif
