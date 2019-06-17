#include "RunHeaderv1.h"

using namespace std;

RunHeaderv1::RunHeaderv1()
  : RunNumber(0)
  , TimeStart(0)
  , TimeStop(0)
  , Bfield(-99999.)
{}

void RunHeaderv1::Reset()
{
  // we don't want to reset the run header for each event
  // so just return here
  return;
}

void RunHeaderv1::identify(ostream& out) const
{
  out << "identify yourself: I am an RunHeaderv1 Object" << endl;
  out << "Run no: " << RunNumber << endl;
  out << "Started at: " << ctime(&TimeStart);
  out << "Ended at:   " << ctime(&TimeStop); 
  out << "B field was " << Bfield << endl;
  return;
}

int RunHeaderv1::isValid() const
{
  return((TimeStart) ? 1:0); // return 1 if TimeStart not zero
}

// this is now cheating, we set the begin run time to the
// time of the event with earliest time until we get it from the data base
void RunHeaderv1::set_TimeStart(time_t start)
{
  if (!TimeStart || TimeStart > start)
    {
      TimeStart = start;
    }
  return;
}

// some more cheating, get the oldest event in fragment
// until we get it from the db
void RunHeaderv1::set_TimeStop(time_t stop)
{
  if (TimeStop < stop)
    {
      TimeStop = stop;
    }
  return;
}

