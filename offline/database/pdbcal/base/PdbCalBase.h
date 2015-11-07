#ifndef PDBCALBASE__H
#define PDBCALBASE__H

#include <TObject.h>

#include <iostream>

class PdbCalBase: public TObject
{
 public: 
  /// ctor
  PdbCalBase(){}
  
  virtual ~PdbCalBase(){}

  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const;

 protected:

  ClassDef(PdbCalBase,1);

};
 
#endif
