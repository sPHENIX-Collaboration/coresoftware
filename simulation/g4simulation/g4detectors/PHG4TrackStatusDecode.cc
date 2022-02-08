#include "PHG4TrackStatusDecode.h"

#include <Geant4/G4TrackStatus.hh>

#include <TSystem.h>

#include <iostream>
#include <map>

using namespace std;
static map<int, string> trackstatus;

string PHG4TrackStatusDecode::GetTrackStatus(const int istatus)
{
  if (trackstatus.empty())
  {
    cout << "filling trackstatus map" << endl;
    trackstatus[fAlive] = "fAlive";
    trackstatus[fStopButAlive] = "fStopButAlive";
    trackstatus[fStopAndKill] = "fStopAndKill";
    trackstatus[fKillTrackAndSecondaries] = "fKillTrackAndSecondaries";
    trackstatus[fSuspend] = "fSuspend";
    trackstatus[fPostponeToNextEvent] = "fPostponeToNextEvent";
  }
  if (trackstatus.find(istatus) == trackstatus.end())
  {
    cout << "could not find status " << istatus << " in trackstatus map" << endl;
    gSystem->Exit(1);
  }
  return trackstatus[istatus];
}
