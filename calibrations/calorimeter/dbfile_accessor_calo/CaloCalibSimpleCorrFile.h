#ifndef CALOCALIB_SIMPLECORRFILE_H
#define CALOCALIB_SIMPLECORRFILE_H

#include <calobase/RawTowerDefs.h>

#include <phool/phool.h>

#include <iostream>
#include <map>
#include <utility>

//
class CaloCalibSimpleCorrFile
{
 public:
  typedef std::map<RawTowerDefs::keytype, float> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  CaloCalibSimpleCorrFile(RawTowerDefs::CalorimeterId caloid = RawTowerDefs::NONE)
    : _caloid(caloid)
  {
  }

  virtual void Open(const std::string &) { PHOOL_VIRTUAL_WARN("Open"); }
  virtual void View() { PHOOL_VIRTUAL_WARN("View"); }
  virtual void ViewReadable() { PHOOL_VIRTUAL_WARN("ViewReadable"); }

  virtual ~CaloCalibSimpleCorrFile() { PHOOL_VIRTUAL_WARN("Destructor"); }

  //void Reset() override;
  //int isValid() const override;
  //void identify(std::ostream &os = std::cout) const override;

  void setCalorimeterID(RawTowerDefs::CalorimeterId caloid) { _caloid = caloid; }
  RawTowerDefs::CalorimeterId getCalorimeterID() { return _caloid; }

  virtual ConstIterator AddCorr(const unsigned int ieta, const unsigned int iphi, float corr) = 0;
  ConstIterator AddCorr(RawTowerDefs::keytype key, float corr)
  {
    _corrs[key] = corr;
    return _corrs.find(key);
  }

  float getCorr(RawTowerDefs::keytype key)
  {
    if (_corrs.find(key) != _corrs.end())
      return _corrs[key];
    else
    {
      std::cout << "calibrations/CaloCalibSimpleCorrFile: "
                << "corr not found for key " << key
                << ", returning -999" << std::endl;
      return -999;
    }
  }

  /*
  const float getCorr(RawTowerDefs::keytype key)
  { 
    if (_corrs.find(key) != _corrs.end()) 
      return _corrs[key];  
    else
      {
	std::cout << "calibrations/CaloCalibSimpleCorrFile: "
		  << "corr not found for key " << key 
		  << ", returning -999" << std::endl;
	return -999; 
      }
  }
  */

  virtual float getCorr(const unsigned int ieta, const unsigned int iphi) = 0;
  //  virtual float getCorr(const unsigned int ieta, const unsigned int iphi) const;

  //  float getCorr(const unsigned int ieta, const unsigned int iphi, const unsigned int il );
  //  const float getCorr(const unsigned int ieta, const unsigned int iphi, const unsigned int il) const;

  unsigned int size() const { return _corrs.size(); }

  ConstRange get_corrs() const
  {
    return make_pair(_corrs.begin(), _corrs.end());
  }

  /*
  Range get_corrs() const override
  {
    return make_pair(_corrs.begin(), _corrs.end());
  }
  */

  Iterator find_corr(RawTowerDefs::keytype id) { return _corrs.find(id); }
  ConstIterator find_corr(RawTowerDefs::keytype id) const { return _corrs.find(id); }

  void clear_corrs() { _corrs.clear(); }

  // protected:
  RawTowerDefs::CalorimeterId _caloid;
  Map _corrs;

  //  ClassDefOverride(CaloCalibSimpleCorrFile, 1);
};

#endif
