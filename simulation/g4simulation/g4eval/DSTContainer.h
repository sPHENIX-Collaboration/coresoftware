// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4EVAL_DSTCONTAINER_H
#define G4EVAL_DSTCONTAINER_H

/*!
 * \file DSTContainer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \modified Christof Roland <cer@mit.edu>
 */

#include <phool/PHObject.h>

//! track evaluation container base class
/*! this is the base class. Does nothing */
class DSTContainer: public PHObject
{
  
  public:
  
  //! constructor
  explicit DSTContainer()
  {}
  
  //! copy constructor
  explicit DSTContainer(const DSTContainer &) = delete;
  
  //! assignment operator
  DSTContainer& operator = ( const DSTContainer& ) = delete;

  ClassDefOverride(DSTContainer,1)
    
};

#endif
