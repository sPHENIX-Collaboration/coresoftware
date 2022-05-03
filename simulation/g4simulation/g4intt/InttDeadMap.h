#ifndef G4INTT_INTTDEADMAP_H
#define G4INTT_INTTDEADMAP_H

#include <g4detectors/PHG4CellDefs.h>

#include <phool/PHObject.h>

#include <iostream>  // for cout, ostream
#include <set>

class InttDeadMap : public PHObject
{
 public:
  typedef std::set<PHG4CellDefs::keytype> Map;

  ~InttDeadMap() override {}

  int isValid() const override;
  void identify(std::ostream &os = std::cout) const override;

  void addDeadChannelIntt(const int layer,
                          const int ladder_phi, const int ladder_z,
                          const int strip_z, const int strip_phi);
  virtual void addDeadChannel(PHG4CellDefs::keytype) { return; }

  virtual bool isDeadChannel(PHG4CellDefs::keytype) const { return false; }
  bool isDeadChannelIntt(const int layer,
                         const int ladder_phi, const int ladder_z,
                         const int strip_z, const int strip_phi) const;

  //! return all towers
  virtual const Map &getDeadChannels(void) const;
  virtual Map &getDeadChannels(void);

  virtual unsigned int size() const { return 0; }

  static PHG4CellDefs::keytype getInttKey(int layer,
                                          int ladder_phi, int ladder_z,
                                          int strip_z, int strip_phi);

  static int getWildCardID() { return s_wildCardID; }

 protected:
  InttDeadMap()
  {
  }

 private:
  static int s_wildCardID;

  ClassDefOverride(InttDeadMap, 1)
};

#endif /* G4INTT_INTTDEADMAP_H */
