#ifndef SYNCOBJECTv2_H
#define SYNCOBJECTv2_H

#include "SyncObjectv1.h"
#include <iostream>

class SyncObjectv2: public SyncObjectv1
{
 public:

  /// ctor
  SyncObjectv2();
  SyncObjectv2(const SyncObject&);

  SyncObjectv2 *clone() const {return new SyncObjectv2(*this);}
  /// dtor
  virtual ~SyncObjectv2() {}

  ///  Clear Event
  void Reset();

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream& os = std::cout) const;

  /// set Segment Number
  void SegmentNumber(const int ival) {segmentnumber = ival;}

 protected:

  /// get Segment Number
  int SegmentNumber() const {return segmentnumber;}

  int segmentnumber; // segment number

 private: // prevent doc++ from showing ClassDef
  ClassDef(SyncObjectv2,1)
};

#endif
