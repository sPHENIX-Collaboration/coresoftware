//-----------------------------------------------------------------------------
//  $Header: /afs/rhic.bnl.gov/phenix/PHENIX_CVS/online_distribution/newbasic/PHTimeStamp.C,v 1.23 2009/08/18 23:03:01 pinkenbu Exp $
//
//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PHTimeStamp
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#include "PHTimeStamp.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

#ifndef HAVE_STRPTIME_PROTOTYPE
extern "C" {
char *strptime(const char *s, const char *format, struct
               tm *tm);
}
#endif

#ifdef WIN32
const phtime_t ticOffset = 35067168000000000UL;
#else
const phtime_t ticOffset = 35067168000000000ULL;
#endif

const phtime_t ticFactor = 10000000;

PHTimeStamp::PHTimeStamp()
{
  setTics(time(0));

  setenv("TZ", "EST5EDT", 1);
}

PHTimeStamp::PHTimeStamp(const int year, const int month, const int day, const int hour, const int minute, const int second, const int fraction)
{
  set(year, month, day, hour, minute, second, fraction);
  setenv("TZ", "EST5EDT", 1);
}

PHTimeStamp::PHTimeStamp(const time_t t)
{
  setTics(t);
  setenv("TZ", "EST5EDT", 1);
}

void PHTimeStamp::set(const int year, const int month, const int day,
                      const int hour, const int minute,
                      const int second, const int fraction)
{
  if (year < 1900)
  {
    setTics(0);
    return;
  }
  tm newTime;
  newTime.tm_year = year - 1900;
  newTime.tm_mon = month - 1;
  newTime.tm_mday = day;
  newTime.tm_hour = hour;
  newTime.tm_min = minute;
  newTime.tm_sec = second;

  // This tells mktime that it's not known whether it's daylight
  // savings time or not.
  newTime.tm_isdst = -1;

  setTics(mktime(&newTime));
  binaryTime += fraction;
}

void PHTimeStamp::set(const char *timeString)
{
  tm newTime;
#ifndef WIN32
  strptime(timeString, "%A %h %d %H:%M:%S %Y", &newTime);
  setTics(mktime(&newTime));
#endif
}

void PHTimeStamp::setToSystemTime()
{
  setTics(time(0));
}

time_t PHTimeStamp::getTics() const
{
  return binaryTimeToTics(binaryTime);
}

void PHTimeStamp::setTics(const time_t tics)
{
  binaryTime = ticsToBinaryTime(tics);
}

void PHTimeStamp::setBinTics(const phtime_t t)
{
  binaryTime = t;
}

phtime_t PHTimeStamp::ticsToBinaryTime(time_t tics) const
{
  return tics * ticFactor + ticOffset;
}

time_t PHTimeStamp::binaryTimeToTics(phtime_t bt) const
{
  return (bt - ticOffset) / ticFactor;
}

int PHTimeStamp::isInRange(const PHTimeStamp &t1, const PHTimeStamp &t2)
{
  return (binaryTime > t1.getBinaryTime() && binaryTime < t2.getBinaryTime());
}

void PHTimeStamp::print()
{
  cout << *this;
}

//
// Operators
//
int PHTimeStamp::operator==(const PHTimeStamp &t) const
{
  return binaryTime == t.getBinaryTime();
}

int PHTimeStamp::operator!=(const PHTimeStamp &t) const
{
  return binaryTime != t.getBinaryTime();
}

int PHTimeStamp::operator>(const PHTimeStamp &t) const
{
  return binaryTime > t.getBinaryTime();
}

int PHTimeStamp::operator<(const PHTimeStamp &t) const
{
  return binaryTime < t.getBinaryTime();
}

int PHTimeStamp::operator>=(const PHTimeStamp &t) const
{
  return binaryTime >= t.getBinaryTime();
}

int PHTimeStamp::operator<=(const PHTimeStamp &t) const
{
  return binaryTime <= t.getBinaryTime();
}

PHTimeStamp &PHTimeStamp::operator=(const PHTimeStamp &t)
{
  binaryTime = t.getBinaryTime();
  return *this;
}

PHTimeStamp PHTimeStamp::operator+=(time_t t)
{
  binaryTime += t * ticFactor;
  return *this;
}

PHTimeStamp PHTimeStamp::operator-=(time_t t)
{
  binaryTime -= t * ticFactor;
  return *this;
}

void PHTimeStamp::print() const
{
  cout << *this << endl;
}
char *PHTimeStamp::formatTimeString() const
{
  // this one gives, for naming purposes, the time string
  // without blanks

  time_t tics = getTics();
  char timeString[25];
  timeString[24] = '\0';
  strncpy(timeString, ctime(&tics), 24);
  char *line = new char[25];

  char *u = strtok(timeString, " ");

  if (u) strcpy(line, u);

  while ((u = strtok(0, " ")))
  {
    strcat(line, "_");
    strcat(line, u);
  }
  return line;
}

//
// Non member operators and functions
//

PHTimeStamp operator+(const PHTimeStamp &t1, time_t t2)
{
  PHTimeStamp newTime = t1;
  newTime += t2;
  return newTime;
}

PHTimeStamp operator-(const PHTimeStamp &t1, time_t t2)
{
  PHTimeStamp newTime = t1;
  newTime -= t2;
  return newTime;
}

time_t operator-(const PHTimeStamp &t1, const PHTimeStamp &t2)
{
  return t1.getTics() - t2.getTics();
}

ostream &operator<<(ostream &s, const PHTimeStamp &t)
{
  time_t tics = t.getTics();
  char timeString[25];
  timeString[24] = '\0';
  strncpy(timeString, ctime(&tics), 24);
  return s << timeString;
}

istream &operator>>(istream &s, PHTimeStamp &t)
{
  char timeString[25];
  s.getline(timeString, 25);
  t.set(timeString);
  return s;
}
