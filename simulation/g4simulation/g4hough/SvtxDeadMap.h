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

  virtual void addDeadChannel(const unsigned int layer, const unsigned int ieta, const unsigned int iphi);
  virtual void addDeadChannel(PHG4CellDefs::keytype key);

  virtual bool isDeadChannel(PHG4CellDefs::keytype key);
  virtual bool isDeadChannel(const unsigned int layer, const unsigned int ieta, const unsigned int iphi);
  //! return all towers
  virtual const Map &getDeadChannels(void) const;
  virtual Map &getDeadChannels(void);

  virtual unsigned int size() const { return 0; }
 protected:
  SvtxDeadMap()
  {
  }

 private:
  ClassDef(SvtxDeadMap, 1)
};

#endif /* SvtxDeadMap_H__ */
