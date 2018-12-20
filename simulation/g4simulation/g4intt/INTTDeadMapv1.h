#ifndef INTTDeadMapv1_H__
#define INTTDeadMapv1_H__

#include "INTTDeadMap.h"

class INTTDeadMapv1 : public INTTDeadMap
{
 public:

  INTTDeadMapv1()
  {
  }
  virtual ~INTTDeadMapv1() {}
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

  ClassDef(INTTDeadMapv1, 1)
};

#endif /* INTTDeadMapv1_H__ */
