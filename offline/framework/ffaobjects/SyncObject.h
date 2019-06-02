// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_SYNCOBJECT_H
#define FFAOBJECTS_SYNCOBJECT_H

#include <phool/PHObject.h>

#include <iostream>

///
class SyncObject: public PHObject
{
 public:

  /// dtor
  virtual ~SyncObject() {}

  /// Clear Sync
  virtual void Reset();

  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const;


  /// isValid returns non zero if object contains valid data
  virtual int isValid() const;

  virtual SyncObject* clone() const;
  virtual SyncObject& operator=(const SyncObject &source);
  virtual int Different(const SyncObject *other) const;

  /// set Event Counter
  virtual void EventCounter(const int /*ival*/) {return;}

  /// set Event Number
  virtual void EventNumber(const int /*ival*/) {return;}


  /// set Segment Number
  virtual void SegmentNumber(const int /*ival*/) {return;}

  /// set Run Number
  virtual void RunNumber(const int /*ival*/) {return;}

 protected:
  /// get Event Number
  virtual int EventNumber() const {return -9999;}
  /// get Event Counter
  virtual int EventCounter() const {return -9999;}
  /// get Run Number
  virtual int RunNumber() const {return -9999;}
  /// get Segment Number
  virtual int SegmentNumber() const {return -9999;}

 private: // prevent doc++ from showing ClassDef
  friend class SyncObjectv1;
  friend class SyncObjectv2;
  friend class Fun4AllDstInputManager;
  friend class DumpSyncObject;
  friend class SegmentSelect;

  ClassDef(SyncObject,1)

};

#endif
