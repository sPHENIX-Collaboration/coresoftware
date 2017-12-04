#ifndef __PHOPERATION_H__
#define __PHOPERATION_H__

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
