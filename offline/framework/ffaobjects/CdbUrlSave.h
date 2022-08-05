// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_CDBURLSAVE_H
#define FFAOBJECTS_CDBURLSAVE_H

#include <phool/PHObject.h>

#include <cstdint>  // for uint64_t
#include <iostream>
#include <string>  // for string
#include <tuple>
#include <vector>  // for vector, vector<>::const_iterator

///
class CdbUrlSave : public PHObject
{
 public:
  /// dtor
  ~CdbUrlSave() override {}

  PHObject *CloneMe() const override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  virtual void AddUrl(const std::string &, const std::string &, const uint64_t) { return; }
  virtual void AddUrl(const std::tuple<std::string, std::string, uint64_t> &) { return; }

  virtual std::vector<std::tuple<std::string, std::string, uint64_t>>::const_iterator begin() const;
  virtual std::vector<std::tuple<std::string, std::string, uint64_t>>::const_iterator end() const;

 private:
  ClassDefOverride(CdbUrlSave, 1)
};

#endif
