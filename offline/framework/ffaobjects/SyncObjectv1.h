#ifndef __SYNCOBJECTv1_H
#define __SYNCOBJECTv1_H

#include "SyncObject.h"
#include <iostream>

class SyncObjectv1: public SyncObject
{
 public:

  /// ctor
  SyncObjectv1();
  SyncObjectv1(const SyncObject&);

  SyncObjectv1 *clone() const {return new SyncObjectv1(*this);}
  /// dtor
  virtual ~SyncObjectv1() {}

  ///  Clear Event
  void Reset();

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream& os = std::cout) const;

  /// isValid returns non zero if object contains valid data
  int isValid() const;

  /// set Event Counter
  void EventCounter(const int ival) {eventcounter = ival;}

  /// set Event Number
  void EventNumber(const int ival) {eventnumber = ival;}

  /// set Run Number
  void RunNumber(const int ival) {runnumber = ival;}

  /// set Segment Number
  void SegmentNumber(const int ival) {segmentnumber = ival;}

 protected:
  /// get Event Counter
  int EventCounter() const {return eventcounter;}
  /// get Event Number
  int EventNumber() const {return eventnumber;}
  /// get Run Number
  int RunNumber() const {return runnumber;}

private:
  int eventcounter;      // running counter
  int eventnumber;  // Event number
  int runnumber;  // Run number
  int segmentnumber; // segment number

  ClassDef(SyncObjectv1,1)
};

#endif
