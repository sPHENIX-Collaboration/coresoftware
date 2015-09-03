#ifndef RAWTOWERCONTAINER_H__
#define RAWTOWERCONTAINER_H__

#include <phool/PHObject.h>
#include <phool/phool.h>
#include <iostream>
#include <map>

class RawTower;

class RawTowerContainer : public PHObject 
{

 public:

  typedef std::map<unsigned int,RawTower *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  RawTowerContainer() {}
  virtual ~RawTowerContainer() {}

  void Reset();
  int isValid() const;
  void identify(std::ostream& os=std::cout) const;
  ConstIterator AddTower(const int ieta, const int iphi, RawTower *twr);
  RawTower *getTower(const int ieta, const int iphi);
  RawTower *getTower(const unsigned int id);
  //! return all towers
  ConstRange getTowers( void ) const;
  Range getTowers( void );

  unsigned int size() const {return _towers.size();}
  void compress(const double emin);
  double getTotalEdep() const;
  unsigned int genkey(const unsigned int ieta, const unsigned int iphi) const;

 protected:
  Map _towers;

  ClassDef(RawTowerContainer,1)
};

#endif /* RAWTOWERCONTAINER_H__ */
