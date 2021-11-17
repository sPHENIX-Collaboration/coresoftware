// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4EVAL_TRACKEVALUATIONCONTAINER_H
#define G4EVAL_TRACKEVALUATIONCONTAINER_H

/*!
 * \file TrackEvaluationContainer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <phool/PHObject.h>

//! track evaluation container base class
/*! this is the base class. Does nothing */
class TrackEvaluationContainer: public PHObject
{
  
  public:
  
  //! constructor
  explicit TrackEvaluationContainer()
  {}
  
  //! copy constructor
  explicit TrackEvaluationContainer(const TrackEvaluationContainer &) = delete;
  
  //! assignment operator
  TrackEvaluationContainer& operator = ( const TrackEvaluationContainer& ) = delete;

  ClassDefOverride(TrackEvaluationContainer,1)
    
};

#endif
