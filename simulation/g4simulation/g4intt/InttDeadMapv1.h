#ifndef G4INTT_INTTDEADMAPV1_H
#define G4INTT_INTTDEADMAPV1_H

#include "InttDeadMap.h"

#include <g4detectors/PHG4CellDefs.h>

#include <iostream>  // for cout, ostream

class InttDeadMapv1 : public InttDeadMap
{
 public:
  InttDeadMapv1()
  {
  }
  ~InttDeadMapv1() override {}

  void Reset() override;
  int isValid() const override;

  void identify(std::ostream &os = std::cout) const override;

  void addDeadChannel(PHG4CellDefs::keytype key) override;

  bool isDeadChannel(PHG4CellDefs::keytype key) const override;
  //! return all towers
  const Map &getDeadChannels(void) const override;
  Map &getDeadChannels(void) override;

  unsigned int size() const override { return m_DeadChannels.size(); }

 private:
  Map m_DeadChannels;

  ClassDefOverride(InttDeadMapv1, 1)
};

#endif /* G4INTT_INTTDEADMAPV1_H */
