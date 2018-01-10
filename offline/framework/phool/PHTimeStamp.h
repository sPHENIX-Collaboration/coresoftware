//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Purpose: PHENIX - wide time stamp class
//
//  Declaration of class PHTimeStamp
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#ifndef PHTIMESTAMP_H__
#define PHTIMESTAMP_H__

#include "PHObject.h"

#include <ctime>

typedef unsigned long long phtime_t;

class PHTimeStamp : public PHObject
{
 public:
  enum
  {
    PHFarFuture = 2147483647
  };

  PHTimeStamp();

  PHTimeStamp(const int, const int, const int, const int, const int, const int, const int = 0);
  PHTimeStamp(const time_t);
  void setBinTics(const phtime_t t);

  virtual ~PHTimeStamp() {}
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

  PHTimeStamp &operator=(const PHTimeStamp &);
  PHTimeStamp operator+=(time_t);
  PHTimeStamp operator-=(time_t);

  char *formatTimeString() const;
  void print() const;

 private:
  phtime_t ticsToBinaryTime(time_t) const;
  time_t binaryTimeToTics(phtime_t) const;

 protected:
  phtime_t binaryTime;
  ClassDef(PHTimeStamp, 1)
};

PHTimeStamp operator+(const PHTimeStamp &, time_t);
PHTimeStamp operator-(const PHTimeStamp &, time_t);
time_t operator-(const PHTimeStamp &, const PHTimeStamp &);

std::ostream &operator<<(std::ostream &, const PHTimeStamp &);
std::istream &operator>>(std::istream &, PHTimeStamp &);

#endif /* __PHTIMESTAMP_H__ */
