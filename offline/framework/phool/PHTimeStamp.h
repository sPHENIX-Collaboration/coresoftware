//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Purpose: PHENIX - wide time stamp class
//
//  Declaration of class PHTimeStamp
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#ifndef PHOOL_PHTIMESTAMP_H
#define PHOOL_PHTIMESTAMP_H

#include "PHObject.h"

#include <ctime>
#include <iosfwd> 

typedef unsigned long long phtime_t;

class PHTimeStamp : public PHObject
{
 public:
  static const unsigned long long PHFarFuture;  // set to ULLONG_MAX

  PHTimeStamp();

  PHTimeStamp(const int, const int, const int, const int, const int, const int, const int = 0);
  PHTimeStamp(const time_t);
  void setBinTics(const phtime_t t);

  ~PHTimeStamp() override {}

 public:
  void set(const int, const int, const int, const int, const int, const int, const int = 0);

  void set(const char *);

  void setToSystemTime();
  void setToFarFuture() { setTics(PHFarFuture); }
  phtime_t getBinaryTime() const { return binaryTime; }
  time_t getTics() const;
  void setTics(const time_t);

  int isInRange(const PHTimeStamp &, const PHTimeStamp &);
  void print();

  int operator==(const PHTimeStamp &) const;
  int operator!=(const PHTimeStamp &) const;
  int operator>(const PHTimeStamp &) const;
  int operator>=(const PHTimeStamp &) const;
  int operator<(const PHTimeStamp &) const;
  int operator<=(const PHTimeStamp &) const;

  PHTimeStamp operator+=(time_t);
  PHTimeStamp operator-=(time_t);

  char *formatTimeString() const;
  void print() const;

 private:
  phtime_t ticsToBinaryTime(time_t) const;
  time_t binaryTimeToTics(phtime_t) const;

 protected:
  phtime_t binaryTime;
  ClassDefOverride(PHTimeStamp, 1)
};

PHTimeStamp operator+(const PHTimeStamp &, time_t);
PHTimeStamp operator-(const PHTimeStamp &, time_t);
time_t operator-(const PHTimeStamp &, const PHTimeStamp &);

std::ostream &operator<<(std::ostream &, const PHTimeStamp &);
std::istream &operator>>(std::istream &, PHTimeStamp &);

#endif /* PHOOL_PHTIMESTAMP_H */
