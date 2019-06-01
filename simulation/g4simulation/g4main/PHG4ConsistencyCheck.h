// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4CONSISTENCYCHECK_H
#define G4MAIN_PHG4CONSISTENCYCHECK_H

#include <fun4all/SubsysReco.h>

#include <string>                // for string

class PHCompositeNode;

class PHG4ConsistencyCheck: public SubsysReco
{
 public:
PHG4ConsistencyCheck( const std::string &name = "CONSISTENCYCHECK" );
 virtual ~PHG4ConsistencyCheck() {}
  //! init
  int InitRun(PHCompositeNode *);

  //! event processing
  int process_event(PHCompositeNode *);

 protected:
  unsigned int errorcnt;

};

#endif
