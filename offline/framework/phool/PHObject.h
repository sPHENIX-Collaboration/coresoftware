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
  virtual ~PHObject() {}
  /// Virtual copy constructor.
  virtual PHObject* CloneMe() const;

#if !defined(__CINT__) || defined(__CLING__)
  virtual PHObject* clone() const final;
  virtual PHObject *Clone(const char *newname = "") const final;
  virtual void 	Copy(TObject &object) const final;
#else
  virtual PHObject* clone() const;
  virtual PHObject *Clone(const char *newname = "") const;
  virtual void 	Copy(TObject &object) const;
#endif

  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const;

  /// Clear Event
  virtual void Reset();

  /// isValid returns non zero if object contains vailid data
  virtual int isValid() const;
  virtual int isValid(const float) const;
  virtual int isValid(const double) const;
  virtual int isValid(const int) const;
  virtual int isValid(const unsigned int) const;

  virtual int isImplemented(const float f) const;
  virtual int isImplemented(const double f) const;
  virtual int isImplemented(const int i) const;
  virtual int isImplemented(const unsigned int i) const;

  virtual int Integrate() const { return 0; }
  virtual int Integrate(PHObject* obj) { return -1; }
  virtual void CopyFrom(const PHObject *obj);

 private:
  ClassDef(PHObject, 0)  // no I/O
};

#endif /* PHOOL_PHOBJECT_H */
