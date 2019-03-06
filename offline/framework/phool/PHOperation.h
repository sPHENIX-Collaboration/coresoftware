#ifndef PHOOL_PHOPERATION_H
#define PHOOL_PHOPERATION_H

//  Declaration of class PHOperation
//  Purpose: abstract strategy base class
//  Author: Matthias Messer

template <class T>
class PHOperation
{
 public:
  PHOperation();
  virtual ~PHOperation();

 public:
  virtual void perform(T*) = 0;
  void
  operator()(T& o)
  {
    perform(&o);
  }
  void
  operator()(T* o)
  {
    perform(o);
  }
};

#endif /* __PHOPERATION_H__ */
