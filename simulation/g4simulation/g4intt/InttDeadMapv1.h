#ifndef G4INTT_INTTDEADMAPV1_H
#define G4INTT_INTTDEADMAPV1_H

#include "InttDeadMap.h"

class InttDeadMapv1 : public InttDeadMap
{
 public:
  InttDeadMapv1()
  {
  }
  virtual ~InttDeadMapv1() {}
  virtual void Reset();
  virtual int isValid() const;

  virtual void identify(std::ostream &os = std::cout) const;
  void addDeadChannel(PHG4CellDefs::keytype key);

  bool isDeadChannel(PHG4CellDefs::keytype key) const;
  //! return all towers
  virtual const Map &getDeadChannels(void) const;
  virtual Map &getDeadChannels(void);

  virtual unsigned int size() const { return m_DeadChannels.size(); }

 private:
  Map m_DeadChannels;

  ClassDef(InttDeadMapv1, 1)
};

#endif /* G4INTT_INTTDEADMAPV1_H */
