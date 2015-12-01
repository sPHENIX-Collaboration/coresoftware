#ifndef PDBCALBASE__H
#define PDBCALBASE__H

#include <phool/PHObject.h>

#include <iostream>

class PdbCalBase: public PHObject
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
