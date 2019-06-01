#ifndef PHG4CONSISTENCYCHECK_H
#define PHG4CONSISTENCYCHECK_H

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
