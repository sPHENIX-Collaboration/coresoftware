#ifndef PHOBJECT_H__
#define PHOBJECT_H__

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
  PHObject();

  /// dtor
  virtual ~PHObject() {}
  /// Virtual copy constructor.
  virtual PHObject* clone() const;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const;

  /// Clear Event
  virtual void Reset();

  /// For integral objects, e.g. integrated luminosity counter, integrate with another object from another run
  virtual void Integrate(PHObject *);

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
  void SplitLevel(const int i) { split = i; }
  int SplitLevel() const { return split; }
  void BufferSize(const int i) { bufSize = i; }
  int BufferSize() const { return bufSize; }
 private:
  int split;    //! not saved, it is set to change the split level for this object
  int bufSize;  //! not saved, it is set to change the buffer size for this object

  ClassDef(PHObject, 0)  // no I/O
};

#endif /* PHOBJECT_H__ */
