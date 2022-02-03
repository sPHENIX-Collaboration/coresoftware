#ifndef PHOOL_PHOBJECT_H
#define PHOOL_PHOBJECT_H

//  Declaration of class PHObject
//  Purpose: Tiny layer between TObject and output objects and
//           enforce some standards
//  Author: Chris Pinkenburg

#include <TObject.h>

#include <iostream>

class PHObject : public TObject
{
 public:
  /// ctor
  PHObject() {}

  /// dtor
  ~PHObject() override {}
  /// Virtual copy constructor.
  virtual PHObject* CloneMe() const;

  virtual PHObject* clone() const final;
  PHObject *Clone(const char *newname = "") const final;
  void 	Copy(TObject &object) const final;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const;

  /// Clear Event
  virtual void Reset();

  /// isValid returns non zero if object contains vailid data
  virtual int isValid() const;

  virtual int Integrate() const { return 0; }
  virtual int Integrate(PHObject* /*obj*/) { return -1; }
  virtual void CopyFrom(const PHObject *obj);

 private:
  ClassDefOverride(PHObject, 0)  // no I/O
};

#endif /* PHOOL_PHOBJECT_H */
