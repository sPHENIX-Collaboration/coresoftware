#include "PHG4TrackStatusDecode.h"

#include <Geant4/G4TrackStatus.hh>

#include <TSystem.h>

#include <iostream>
#include <map>

static std::map<int, std::string> trackstatus;

std::string PHG4TrackStatusDecode::GetTrackStatus(const int istatus)
{
  if (trackstatus.empty())
  {
    std::cout << "filling trackstatus map" << std::endl;
    trackstatus[fAlive] = "fAlive";
    trackstatus[fStopButAlive] = "fStopButAlive";
    trackstatus[fStopAndKill] = "fStopAndKill";
    trackstatus[fKillTrackAndSecondaries] = "fKillTrackAndSecondaries";
    trackstatus[fSuspend] = "fSuspend";
    trackstatus[fPostponeToNextEvent] = "fPostponeToNextEvent";
  }
  if (trackstatus.find(istatus) == trackstatus.end())
  {
    std::cout << "could not find status " << istatus << " in trackstatus map" << std::endl;
    gSystem->Exit(1);
  }
  return trackstatus[istatus];
}
