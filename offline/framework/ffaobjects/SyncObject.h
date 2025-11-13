// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_SYNCOBJECT_H
#define FFAOBJECTS_SYNCOBJECT_H

#include <phool/PHObject.h>

#include <iostream>
#include <limits>

///
class SyncObject : public PHObject
{
 public:
  /// ctor - daughter class copy ctor needs this
  SyncObject() = default;
  /// copy ctor daughter class copy ctor needs also this
  SyncObject(const SyncObject& source) = default;
  /// dtor
  ~SyncObject() override = default;
  /// Clear Sync
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  PHObject* CloneMe() const override;
  SyncObject& operator=(const SyncObject& source);
  virtual unsigned int Different(const SyncObject* other) const;

  /// set Event Counter
  virtual void EventCounter(const int /*ival*/) { return; }

  /// set Event Number
  virtual void EventNumber(const int /*ival*/) { return; }

  /// set Segment Number
  virtual void SegmentNumber(const int /*ival*/) { return; }

  /// set Run Number
  virtual void RunNumber(const int /*ival*/) { return; }

  /// get Event Number
  virtual int EventNumber() const { return std::numeric_limits<int>::min(); }

protected:
  /// get Event Counter
  virtual int EventCounter() const { return std::numeric_limits<int>::min(); }
  /// get Run Number
  virtual int RunNumber() const { return std::numeric_limits<int>::min(); }
  /// get Segment Number
  virtual int SegmentNumber() const { return std::numeric_limits<int>::min(); }

 private:  // prevent doc++ from showing ClassDefOverride
  friend class SyncObjectv1;
  friend class Fun4AllDstInputManager;
  friend class Fun4AllDstPileupInputManager;
  friend class DumpSyncObject;
  friend class SegmentSelect;

  ClassDefOverride(SyncObject, 1)
};

#endif
