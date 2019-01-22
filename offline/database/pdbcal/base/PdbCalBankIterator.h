#ifndef PDBCAL_BASE_PDBCALBANKITERATOR_H
#define PDBCAL_BASE_PDBCALBANKITERATOR_H

class PdbBankID;
class PdbCalBank;

#include <phool/PHTimeStamp.h>

#include <iostream>

class PdbCalBankIterator 
{
public:
  virtual ~PdbCalBankIterator() {}

  virtual bool init(const std::string &fulldbname, const PdbBankID& bankid) = 0;
  
  virtual bool isValid() const = 0;

  virtual void print(std::ostream& os = std::cout) const = 0;

  virtual void setBankID(const PdbBankID& id) = 0;

  virtual void setEndValTimeLimits
  (const PHTimeStamp& min = PHTimeStamp(0),
   const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture)) = 0;
  
  virtual void setInsertTimeLimits
  (const PHTimeStamp& min = PHTimeStamp(0),
   const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture)) = 0;
  
  virtual void setStartValTimeLimits
  (const PHTimeStamp& min = PHTimeStamp(0),
   const PHTimeStamp& max = PHTimeStamp(PHTimeStamp::PHFarFuture)) = 0;

  virtual PdbCalBank* next() = 0;
};

#endif
