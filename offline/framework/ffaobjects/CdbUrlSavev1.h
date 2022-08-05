// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_CDBURLSAVEV1_H
#define FFAOBJECTS_CDBURLSAVEV1_H

#include "CdbUrlSave.h"

#include <cstdint>  // for uint64_t
#include <iostream>
#include <string>  // for string
#include <tuple>
#include <vector>  // for vector<>::const_iterator, vector

///
class CdbUrlSavev1 : public CdbUrlSave
{
 public:
  /// dtor
  ~CdbUrlSavev1() override {}

  PHObject *CloneMe() const override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  void AddUrl(const std::string &domain, const std::string &url, const uint64_t timestamp) override;
  void AddUrl(const std::tuple<std::string, std::string, uint64_t> &tup) override;

  std::vector<std::tuple<std::string, std::string, uint64_t>>::const_iterator begin() const override;
  std::vector<std::tuple<std::string, std::string, uint64_t>>::const_iterator end() const override;

 private:
  std::vector<std::tuple<std::string, std::string, uint64_t>> m_CdbUrlVector;

  ClassDefOverride(CdbUrlSavev1, 1)
};

#endif
