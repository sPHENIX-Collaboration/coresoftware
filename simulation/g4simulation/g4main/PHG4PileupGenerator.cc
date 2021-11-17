#include "PHG4PileupGenerator.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <gsl/gsl_randist.h>

#include <iostream>  // for operator<<, endl, basic_ostream

using namespace std;

PHG4PileupGenerator::PHG4PileupGenerator(const string &name)
  : PHG4ParticleGeneratorBase(name)
{
  return;
}

PHG4PileupGenerator::~PHG4PileupGenerator()
{
  delete _generator;
}

int PHG4PileupGenerator::Init(PHCompositeNode *topNode)
{
  if (!_generator) return Fun4AllReturnCodes::EVENT_OK;
  _generator->Init(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4PileupGenerator::InitRun(PHCompositeNode *topNode)
{
  //cout << "generator ptr: " << _generator << endl;

  if (!_generator) return Fun4AllReturnCodes::EVENT_OK;

  _generator->InitRun(topNode);

  _ave_coll_per_crossing = _collision_rate * _time_between_crossings * 1000.0 * 1e-9;

  _min_crossing = _min_integration_time / _time_between_crossings;
  _max_crossing = _max_integration_time / _time_between_crossings;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4PileupGenerator::process_event(PHCompositeNode *topNode)
{
  if (!_generator) return Fun4AllReturnCodes::EVENT_OK;

  // toss multiple crossings all the way back
  for (int icrossing = _min_crossing; icrossing <= _max_crossing; ++icrossing)
  {
    double crossing_time = _time_between_crossings * icrossing;

    int ncollisions = gsl_ran_poisson(RandomGenerator(), _ave_coll_per_crossing);
    if (icrossing == 0) --ncollisions;

    for (int icollision = 0; icollision < ncollisions; ++icollision)
    {
      _generator->set_t0(crossing_time);
      _generator->process_event(topNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4PileupGenerator::Reset(PHCompositeNode *topNode)
{
  if (!_generator) return Fun4AllReturnCodes::EVENT_OK;
  _generator->Reset(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4PileupGenerator::ResetEvent(PHCompositeNode *topNode)
{
  if (!_generator) return Fun4AllReturnCodes::EVENT_OK;
  _generator->ResetEvent(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4PileupGenerator::EndRun(const int runnumber)
{
  if (!_generator) return Fun4AllReturnCodes::EVENT_OK;
  _generator->EndRun(runnumber);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4PileupGenerator::End(PHCompositeNode *topNode)
{
  if (!_generator) return Fun4AllReturnCodes::EVENT_OK;
  _generator->End(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}
