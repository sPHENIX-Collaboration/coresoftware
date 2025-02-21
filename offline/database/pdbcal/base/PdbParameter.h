//  Declaration of class PdbParameter
//  Purpose: single parameter storage class
//  Author: federica

#ifndef PDBCAL_BASE_PDBPARAMETER_H
#define PDBCAL_BASE_PDBPARAMETER_H

#include "PdbCalChan.h"

#include <limits>
#include <string>

class PdbParameter : public PdbCalChan
{
 public:
  PdbParameter(const double, const std::string &name);
  ~PdbParameter() override = default;

  double getParameter() const { return thePar; }
  const std::string &getName() const { return theName; }

  void setParameter(const double val) { thePar = val; }
  void setName(const std::string &name) { theName = name; }

  void print() const override;

 protected:
  PdbParameter() = default;  // this ctor should not be called

 private:
  double thePar{std::numeric_limits<double>::quiet_NaN()};
  std::string theName;

  ClassDefOverride(PdbParameter, 1);
};

#endif /* PDBCAL_BASE_PDBPARAMETER_H */
