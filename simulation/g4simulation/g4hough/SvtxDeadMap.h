#ifndef SvtxDeadMap_H__
#define SvtxDeadMap_H__

#include <g4detectors/PHG4CellDefs.h>

#include <phool/PHObject.h>
#include <set>

class SvtxDeadMap : public PHObject
{
 public:
  typedef std::set<PHG4CellDefs::keytype> Map;

  virtual ~SvtxDeadMap() {}
  virtual void Reset();
  virtual int isValid() const;

  virtual void identify(std::ostream &os = std::cout) const;

  void addDeadChannel(const unsigned int layer, const unsigned int ieta, const unsigned int iphi);
  void addDeadChannelINTT(const unsigned int layer,
                                  const unsigned int ladder_z, const unsigned int ladder_phi,
                                  const unsigned int strip_z, const unsigned int strip_phi);
  virtual void addDeadChannel(PHG4CellDefs::keytype key);

  virtual bool isDeadChannel(PHG4CellDefs::keytype key);
  bool isDeadChannel(const unsigned int layer, const unsigned int ieta, const unsigned int iphi);
  bool isDeadChannel(const unsigned int layer,
                             const unsigned int ladder_z, const unsigned int ladder_phi,
                             const unsigned int strip_z, const unsigned int strip_phi) ;

  //! return all towers
  virtual const Map &getDeadChannels(void) const;
  virtual Map &getDeadChannels(void);

  virtual unsigned int size() const { return 0; }

  static PHG4CellDefs::keytype getINTTKey(const unsigned int layer,
                                          const unsigned int ladder_z, const unsigned int ladder_phi,
                                          const unsigned int strip_z, const unsigned int strip_phi);

  static PHG4CellDefs::keytype s_wildCardID;

 protected:
  SvtxDeadMap()
  {
  }

 private:
  ClassDef(SvtxDeadMap, 1)
};

#endif /* SvtxDeadMap_H__ */
