#ifndef PHOOL_PHNODEOPERATION_H
#define PHOOL_PHNODEOPERATION_H

//  Declaration of class PHNodeOperation
//  Purpose: abstract strategy base class which operates on PHNodes
//  Author: Matthias Messer

class PHNode;

class PHNodeOperation
{
 public:
  PHNodeOperation()
    : verbosity(0)
  {
  }
  virtual ~PHNodeOperation() {}
  void
  operator()(PHNode& o)
  {
    perform(&o);
  }
  void
  operator()(PHNode* o)
  {
    perform(o);
  }

  virtual void Verbosity(const int i) { verbosity = i; }
  virtual int Verbosity() const { return verbosity; }

 protected:
  virtual void perform(PHNode*) = 0;
  int verbosity;
};

#endif
