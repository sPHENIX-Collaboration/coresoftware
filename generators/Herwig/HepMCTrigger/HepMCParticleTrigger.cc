#include "HepMCJetTrigger.h"

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

//____________________________________________________________________________..
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
HepMCJetTrigger::HepMCJetTrigger(float trigger_thresh, int n_incom, bool up_lim, const std::string& name)
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
}

//____________________________________________________________________________..
int HepMCJetTrigger::process_event(PHCompositeNode* topNode)
{
  // std::cout << "HepMCJetTrigger::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  n_evts++;
  if (this->set_event_limit == true)
  {  // needed to keep all HepMC output at the same number of events
    if (n_good >= this->goal_event_number)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
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
    bool const good_event = isGoodEvent(ev);
    if (good_event)
    {
      n_good++;
    }
    if (!good_event)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
void HepMCParticleTrigger::AddParticles(const std::string &particles)
{
   std::vector<int> addedParts = convertToInt()
  _theParticles.insert(
}
bool HepMCParticleTrigger::isGoodEvent(HepMC::GenEvent* e1)
{
  // this is really just the call to actually evaluate and return the filter
  if (this->threshold == 0)
  {
    return true;
  }
  std::vector<int> n_trigger_particles = jetsAboveThreshold(jets);
  for(auto ntp:n_trigger_particles)
  {
    if(ntp <=0 ) return false;
  }
  return true;
}

std::vector<int> HepMCParticleTrigger::getParticles(HepMC::GenEvent* e1)
{
  for (HepMC::GenEvent::particle_const_iterator iter = e1->particles_begin(); iter != e1->particles_end(); ++iter)
  {
    if (m_doStableParticleOnly && ((*iter)->end_vertex() && (*iter)->status() != 1)) continue;
    else{
      auto p = (*iter)->momentum(); 
      float px = p.px();
      float py = p.py();
      float pz = p.pz();
      float p = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
      float pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
      int pid = (*iter)->pid();
      double eta = p.eta();
      if((_doEtaHighCut || _doBothEtaCut ) && eta > _theEtaHigh) continue;
      if((_doEtaLowCut || _doBothEtaCut ) && eta < _theEtaLow) continue;
      if((_doAbsEtaHighCut || _doBothAbsEtaCut ) && std::abs(eta) > _theAbsEtaHigh) continue;
      if((_doAbsEtaLowCut || _doBothAbsEtaCut ) && std::abs(eta) < _theAbsEtaLow) continue;
      if((_doPtHighCut || _doBothPtCut ) && pt > _thePtHigh) continue;
      if((_doPLowCut || _doBothPCut ) && pt < _thePtLow) continue;
      if((_doPHighCut || _doBothPCut ) && p > _thePHigh) continue;
      if((_doPLowCut || _doBothPCut ) && p < _thePLow) continue;
      if((_doPzHighCut || _doBothPzCut ) && pz > _thePzHigh) continue;
      if((_doPzLowCut || _doBothPzCut ) && pz < _thePzLow) continue;

   } 
  return output;
}
int HepMCParticle::particlesAboveThreshold(const std::vector<fastjet::PseudoJet>& jets)
{
  // search through for the number of identified jets above the threshold
  int n_good_jets = 0;
  for (const auto& j : jets)
  {
    float const pt = j.pt();
    if (pt > this->threshold)
    {
      n_good_jets++;
    }
  }
  return n_good_jets;
}
