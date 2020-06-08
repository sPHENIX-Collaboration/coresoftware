#ifndef G4EVAL_EVENTCOUNTER_HP_H
#define G4EVAL_EVENTCOUNTER_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHTimer.h>

#include <memory>

class EventCounter_hp : public SubsysReco
{
  public:

  /// constructor
  EventCounter_hp( const std::string& = "EVENTCOUNTER_HP", unsigned int granularity = 100 );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*)  override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  private:

  // event counter
  unsigned int _ievent = 0;
  unsigned int _granularity = 100;
  std::unique_ptr<PHTimer> _timer;

};

#endif
