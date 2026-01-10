#include "HepMCParticleTrigger.h"

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <fastjet/JetDefinition.hh>

#include <HepMC/GenEvent.h>

#include <fastjet/PseudoJet.hh>
#include <string>
#include <vector>
#include <map> 
//____________________________________________________________________________..
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
HepMCParticleTrigger::HepMCParticleTrigger(float trigger_thresh, int n_incom, bool up_lim, const std::string& name)
  : SubsysReco(name)
  , threshold(trigger_thresh)
  , goal_event_number(n_incom)
  , set_event_limit(up_lim)
  , _theEtaHigh(-999.9)
  , _theEtaLow(-999.9)
  , _thePtHigh(999.9)
  , _thePtLow(-999.9)
  , _thePHigh(999.9)
  , _thePLow(-999.9)
  , _thePzHigh(999.9)
  , _thePzLow(-999.9)
  ,

    _doEtaHighCut(false)
  , _doEtaLowCut(false)
  , _doBothEtaCut(false)
  ,

  _doAbsEtaHighCut(false)
  , _doAbsEtaLowCut(false)
  , _doBothAbsEtaCut(false)
  ,

  _doPtHighCut(false)
  , _doPtLowCut(false)
  , _doBothPtCut(false)
  ,

  _doPHighCut(false)
  , _doPLowCut(false)
  , _doBothPCut(false)
  ,

  _doPzHighCut(false)
  , _doPzLowCut(false)
  , _doBothPzCut(false) 
{
	if(threshold != 0 )
	{
		_doPtLowCut=true;
		_thePtLow=threshold;
	}
}

//____________________________________________________________________________..
int HepMCParticleTrigger::process_event(PHCompositeNode* topNode)
{
  // std::cout << "HepMCParticleTrigger::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  n_evts++;
  if (this->set_event_limit == true)
  {  // needed to keep all HepMC output at the same number of events
    if (n_good >= this->goal_event_number)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  bool good_event;
  PHHepMCGenEventMap* phg = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!phg)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  for (PHHepMCGenEventMap::ConstIter eventIter = phg->begin(); eventIter != phg->end(); ++eventIter)
  {
    PHHepMCGenEvent* hepev = eventIter->second;
    if (!hepev)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    HepMC::GenEvent* ev = hepev->getEvent();
    if (!ev)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    good_event = isGoodEvent(ev);
    if (!good_event)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
   if (good_event)
   {
     n_good++;
   }
  return Fun4AllReturnCodes::EVENT_OK;
}
void HepMCParticleTrigger::AddParticle(int particlePid)
{
	_theParticles.push_back(particlePid);
	return;
}
void HepMCParticleTrigger::AddParticles(std::vector<int> particles)
{
	for(auto p:particles) _theParticles.push_back(p);
	return;
}

void HepMCParticleTrigger::SetPtHigh(double pt) 
{
	_thePtHigh=pt;
	_doPtHighCut=true;
	if(_doPtLowCut) _doBothPtCut=true;
	return;
}
void HepMCParticleTrigger::SetPtLow(double pt) 
{
	_thePtLow=pt;
	_doPtLowCut=true;
	if(_doPtHighCut) _doBothPtCut=true;
	return;
}
void HepMCParticleTrigger::SetPtHighLow(double ptHigh, double ptLow) 
{
	_thePtHigh=ptHigh;
	_doPtHighCut=true;
	_thePtLow=ptLow;
	_doPtLowCut=true;
	_doBothPtCut=true;
	return;
}
void HepMCParticleTrigger::SetPHigh(double pt) 
{
	_thePHigh=pt;
	_doPHighCut=true;
	if(_doPLowCut) _doBothPCut=true;
	return;
}
void HepMCParticleTrigger::SetPLow(double pt) 
{
	_thePLow=pt;
	_doPLowCut=true;
	if(_doPHighCut) _doBothPCut=true;
	return;
}
void HepMCParticleTrigger::SetPHighLow(double ptHigh, double ptLow) 
{
	_thePHigh=ptHigh;
	_doPHighCut=true;
	_thePLow=ptLow;
	_doPLowCut=true;
	_doBothPCut=true;
	return;
}
void HepMCParticleTrigger::SetPzHigh(double pt) 
{
	_thePzHigh=pt;
	_doPzHighCut=true;
	if(_doPzLowCut) _doBothPzCut=true;
	return;
}
void HepMCParticleTrigger::SetPzLow(double pt) 
{
	_thePzLow=pt;
	_doPzLowCut=true;
	if(_doPzHighCut) _doBothPzCut=true;
	return;
}
void HepMCParticleTrigger::SetPzHighLow(double ptHigh, double ptLow) 
{
	_thePzHigh=ptHigh;
	_doPzHighCut=true;
	_thePzLow=ptLow;
	_doPzLowCut=true;
	_doBothPzCut=true;
	return;
}
void HepMCParticleTrigger::SetEtaHigh(double pt) 
{
	_theEtaHigh=pt;
	_doEtaHighCut=true;
	if(_doEtaLowCut) _doBothEtaCut=true;
	return;
}
void HepMCParticleTrigger::SetEtaLow(double pt) 
{
	_theEtaLow=pt;
	_doEtaLowCut=true;
	if(_doEtaHighCut) _doBothEtaCut=true;
	return;
}
void HepMCParticleTrigger::SetEtaHighLow(double ptHigh, double ptLow) 
{
	_theEtaHigh=ptHigh;
	_doEtaHighCut=true;
	_theEtaLow=ptLow;
	_doEtaLowCut=true;
	_doBothEtaCut=true;
	return;
}
void HepMCParticleTrigger::SetAbsEtaHigh(double pt) 
{
	_theEtaHigh=pt;
	_doAbsEtaHighCut=true;
	if(_doAbsEtaLowCut) _doBothAbsEtaCut=true;
	return;
}
void HepMCParticleTrigger::SetAbsEtaLow(double pt) 
{
	_theEtaLow=pt;
	_doAbsEtaLowCut=true;
	if(_doAbsEtaHighCut) _doBothAbsEtaCut=true;
	return;
}
void HepMCParticleTrigger::SetAbsEtaHighLow(double ptHigh, double ptLow) 
{
	_theEtaHigh=ptHigh;
	_doAbsEtaHighCut=true;
	_theEtaLow=ptLow;
	_doAbsEtaLowCut=true;
	_doBothAbsEtaCut=true;
	return;
}
bool HepMCParticleTrigger::isGoodEvent(HepMC::GenEvent* e1)
{
  // this is really just the call to actually evaluate and return the filter
  /*if (this->threshold == 0)
  {
    return true;
  }*/
  std::vector<int> n_trigger_particles = getParticles(e1);
  for(auto ntp:n_trigger_particles)
  {
    if(ntp <=0 ) return false; //make sure all particles have at least 1 
  }
  return true;
}

std::vector<int> HepMCParticleTrigger::getParticles(HepMC::GenEvent* e1)
{  
 std::vector<int> n_trigger {};
 std::map<int, int> particle_types;
 for (HepMC::GenEvent::particle_const_iterator iter = e1->particles_begin(); iter != e1->particles_end(); ++iter)
  {
    if (m_doStableParticleOnly && ((*iter)->end_vertex() || (*iter)->status() != 1)) continue;
    else{
      auto p = (*iter)->momentum(); 
      float px = p.px();
      float py = p.py();
      float pz = p.pz();
      float p_M = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
      float pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
      int pid = (*iter)->pdg_id();
      double eta = p.eta();
      if((_doEtaHighCut || _doBothEtaCut ) && eta > _theEtaHigh) continue;
      if((_doEtaLowCut || _doBothEtaCut ) && eta < _theEtaLow) continue;
      if((_doAbsEtaHighCut || _doBothAbsEtaCut ) && std::abs(eta) > _theEtaHigh) continue;
      if((_doAbsEtaLowCut || _doBothAbsEtaCut ) && std::abs(eta) < _theEtaLow) continue;
      if((_doPtHighCut || _doBothPtCut ) && pt > _thePtHigh) continue;
      if((_doPtLowCut || _doBothPtCut ) && pt < _thePtLow) continue;
      if((_doPHighCut || _doBothPCut ) && p_M > _thePHigh) continue;
      if((_doPLowCut || _doBothPCut ) && p_M < _thePLow) continue;
      if((_doPzHighCut || _doBothPzCut ) && pz > _thePzHigh) continue;
      if((_doPzLowCut || _doBothPzCut ) && pz < _thePzLow) continue;
      if(particle_types.find(pid) != particle_types.end()) particle_types[pid]++;
      else particle_types[pid]=1;
     }
  }
  for(auto p:_theParticles)
  {
    n_trigger.push_back(particleAboveThreshold(particle_types, p)); //make sure we have at least one of each required particle
  } 
  return n_trigger;
}
int HepMCParticleTrigger::particleAboveThreshold(const std::map<int, int>& n_particles, int trigger_particle )
{
 // search through for the number of identified trigger particles passing cuts
  auto it = n_particles.find(std::abs(trigger_particle));
  if( it!= n_particles.end()) return it->second;
  else return 0;
}
