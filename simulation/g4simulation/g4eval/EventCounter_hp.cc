#include "EventCounter_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimer.h>

//_____________________________________________________________________
EventCounter_hp::EventCounter_hp( const std::string& name, unsigned int granularity ):
  SubsysReco( name),
  _granularity( granularity )
{ std::cout << "EventCounter_hp::EventCounter_hp." << std::endl; }

//_____________________________________________________________________
int EventCounter_hp::Init(PHCompositeNode* topNode )
{

  std::cout << "EventCounter_hp::Init." << std::endl;

  // initialize timer
  _timer.reset( new PHTimer("_eventCounter_hp_timer") );
  _timer->restart();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int EventCounter_hp::InitRun(PHCompositeNode* )
{
  std::cout << "EventCounter_hp::InitRun." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int EventCounter_hp::process_event(PHCompositeNode* topNode)
{
  // print event number
  if( _granularity > 0 && (_ievent % _granularity) == 0 )
  { std::cout << "EventCounter_hp::process_event - Event = " << _ievent << std::endl; }
  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int EventCounter_hp::End(PHCompositeNode* )
{
  std::cout << "EventCounter_hp::End." << std::endl;

  // print timer information
  _timer->stop();
  std::cout
    << "EventCounter_hp::End -"
    << " time per event:" << _timer->get_accumulated_time()/(1000.*_ievent) << " sec"
    << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
