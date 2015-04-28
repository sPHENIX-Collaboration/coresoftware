#ifndef __PHNODEOPERATION_H__
#define __PHNODEOPERATION_H__

//  Declaration of class PHNodeOperation
//  Purpose: abstract strategy base class which operates on PHNodes
//  Author: Matthias Messer

class PHNode;

class PHNodeOperation 
{ 
public: 
  PHNodeOperation(){} 
  virtual ~PHNodeOperation(){}
  
  void 
  operator () (PHNode& o) 
  { 
    perform(&o); 
  }
  void 
  operator () (PHNode* o) 
  { 
    perform(o); 
  }
  
protected: 
  virtual void perform(PHNode*) = 0;
}; 

#endif /* __PHNODEOPERATION_H__ */
