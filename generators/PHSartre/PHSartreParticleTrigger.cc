#include "PHSartreParticleTrigger.h"

#include <sartre/Event.h>

#include <TLorentzVector.h>  // for TLorentzVector

#include <cmath>            // for abs
#include <iostream>          // for operator<<, endl, basic_ostream, basic_o...

//__________________________________________________________
PHSartreParticleTrigger::PHSartreParticleTrigger(const std::string &name)
  : PHSartreGenTrigger(name)
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

PHSartreParticleTrigger::~PHSartreParticleTrigger()
{
  if (Verbosity() > 0) PrintConfig();
}

bool PHSartreParticleTrigger::Apply(Event *event)
{
  if (Verbosity() > 2)
  {
    cout << "PHSartreParticleTrigger::Apply - event size: "
         << event->particles.size() << endl;
  }

  // Loop over all particles in the event
  for (unsigned int i = 0; i < event->particles.size(); ++i)
  {
    //cout << i << " " << event->particles[i].pdgId << endl;

    if ((i < 7) && (i != 4)) continue;  // only the VM, daughters OR the nuclear remnants

    const Particle &particle = event->particles[i];

    // loop over all the trigger particle criteria
    for (int j = 0; j < int(_theParticles.size()); j++)
    {
      if ((particle.pdgId == _theParticles[j]) &&
          ((particle.status == 1) || (i == 4)))
      {  //only stable particles OR the vector meson

        if (_doBothEtaCut && (particle.p.Eta() < _theEtaLow ||
                              particle.p.Eta() > _theEtaHigh)) continue;
        if (_doEtaLowCut && particle.p.Eta() < _theEtaLow) continue;
        if (_doEtaHighCut && particle.p.Eta() > _theEtaHigh) continue;

        if (_doBothAbsEtaCut && (abs(particle.p.Eta()) < _theEtaLow ||
                                 abs(particle.p.Eta()) > _theEtaHigh)) continue;
        if (_doAbsEtaLowCut && abs(particle.p.Eta()) < _theEtaLow) continue;
        if (_doAbsEtaHighCut && abs(particle.p.Eta()) > _theEtaHigh) continue;

        if (_doBothPtCut && (particle.p.Pt() < _thePtLow ||
                             particle.p.Pt() > _thePtHigh)) continue;
        if (_doPtHighCut && particle.p.Pt() > _thePtHigh) continue;
        if (_doPtLowCut && particle.p.Pt() < _thePtLow) continue;

        if (_doBothPCut && (particle.p.P() < _thePLow ||
                            particle.p.P() > _thePHigh)) continue;
        if (_doPHighCut && particle.p.P() > _thePHigh) continue;
        if (_doPLowCut && particle.p.P() < _thePLow) continue;

        if (_doBothPzCut && (particle.p.Pz() < _thePzLow ||
                             particle.p.Pz() > _thePzHigh)) continue;
        if (_doPzHighCut && particle.p.Pz() > _thePzHigh) continue;
        if (_doPzLowCut && particle.p.Pz() < _thePzLow) continue;

        //If we made it here, success!
        return true;

      }  //if _theParticles
    }    //_theParticles for loop
  }

  return false;
}

void PHSartreParticleTrigger::AddParticles(const std::string &particles)
{
  std::vector<int> addedParts = convertToInts(particles);
  _theParticles.insert(_theParticles.end(), addedParts.begin(), addedParts.end());
}

void PHSartreParticleTrigger::AddParticles(int particle)
{
  _theParticles.push_back(particle);
}

void PHSartreParticleTrigger::AddParticles(std::vector<int> particles)
{
  _theParticles.insert(_theParticles.end(), particles.begin(), particles.end());
}

void PHSartreParticleTrigger::SetPtHigh(double pt)
{
  _thePtHigh = pt;
  if (_doPtLowCut)
    _doBothPtCut = true;
  else
    _doPtHighCut = true;
}

void PHSartreParticleTrigger::SetPtLow(double pt)
{
  _thePtLow = pt;
  if (_doPtHighCut)
    _doBothPtCut = true;
  else
    _doPtLowCut = true;
}

void PHSartreParticleTrigger::SetPtHighLow(double ptHigh, double ptLow)
{
  if (ptHigh < ptLow)
  {
    _thePtHigh = ptLow;
    _thePtLow = ptHigh;
  }
  else
  {
    _thePtHigh = ptHigh;
    _thePtLow = ptLow;
  }
  _doBothPtCut = true;
  _doPtLowCut = false;
  _doPtHighCut = false;
}

void PHSartreParticleTrigger::SetPHigh(double p)
{
  _thePHigh = p;
  if (_doPLowCut)
  {
    _doBothPCut = true;
    _doPLowCut = false;
  }
  else
  {
    _doPHighCut = true;
  }
}

void PHSartreParticleTrigger::SetPLow(double p)
{
  _thePLow = p;
  if (_doPHighCut)
  {
    _doBothPCut = true;
    _doPHighCut = false;
  }
  else
  {
    _doPLowCut = true;
  }
}

void PHSartreParticleTrigger::SetPHighLow(double pHigh, double pLow)
{
  if (pHigh < pLow)
  {
    _thePHigh = pLow;
    _thePLow = pHigh;
  }
  else
  {
    _thePHigh = pHigh;
    _thePLow = pLow;
  }
  _doBothPCut = true;
  _doPLowCut = false;
  _doPHighCut = false;
}

void PHSartreParticleTrigger::SetEtaHigh(double eta)
{
  _theEtaHigh = eta;
  if (_doEtaLowCut)
  {
    _doBothEtaCut = true;
    _doEtaLowCut = false;
  }
  else
  {
    _doEtaHighCut = true;
  }
}

void PHSartreParticleTrigger::SetEtaLow(double eta)
{
  _theEtaLow = eta;
  if (_doEtaHighCut)
  {
    _doBothEtaCut = true;
    _doEtaHighCut = false;
  }
  else
  {
    _doEtaLowCut = true;
  }
}

void PHSartreParticleTrigger::SetEtaHighLow(double etaHigh, double etaLow)
{
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;
  _doBothEtaCut = true;
  _doEtaHighCut = false;
  _doEtaLowCut = false;
}

void PHSartreParticleTrigger::SetAbsEtaHigh(double eta)
{
  _theEtaHigh = eta;
  if (_doAbsEtaLowCut)
  {
    _doBothAbsEtaCut = true;
    _doAbsEtaLowCut = false;
  }
  else
  {
    _doAbsEtaHighCut = true;
  }
}

void PHSartreParticleTrigger::SetAbsEtaLow(double eta)
{
  _theEtaLow = eta;
  if (_doAbsEtaHighCut)
  {
    _doBothAbsEtaCut = true;
    _doAbsEtaHighCut = false;
  }
  else
  {
    _doAbsEtaLowCut = true;
  }
}

void PHSartreParticleTrigger::SetAbsEtaHighLow(double etaHigh, double etaLow)
{
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;
  _doBothAbsEtaCut = true;
  _doAbsEtaLowCut = false;
  _doAbsEtaHighCut = false;
}

void PHSartreParticleTrigger::SetPzHigh(double pz)
{
  _thePzHigh = pz;
  if (_doPzLowCut)
  {
    _doBothPzCut = true;
    _doPzLowCut = false;
  }
  else
  {
    _doPzHighCut = true;
  }
}

void PHSartreParticleTrigger::SetPzLow(double pz)
{
  _thePzLow = pz;
  if (_doPzHighCut)
  {
    _doBothPzCut = true;
    _doPzHighCut = false;
  }
  else
  {
    _doPzLowCut = true;
  }
}

void PHSartreParticleTrigger::SetPzHighLow(double pzHigh, double pzLow)
{
  if (pzHigh < pzLow)
  {
    _thePzHigh = pzLow;
    _thePzLow = pzHigh;
  }
  else
  {
    _thePzHigh = pzHigh;
    _thePzLow = pzLow;
  }
  _doBothPzCut = true;
  _doPzLowCut = false;
  _doPzHighCut = false;
}

void PHSartreParticleTrigger::PrintConfig()
{
  cout << "---------------- PHSartreParticleTrigger::PrintConfig --------------------" << endl;
  cout << "   Particles: ";
  for (int i = 0; i < int(_theParticles.size()); i++) cout << _theParticles[i] << "  ";
  cout << endl;

  if (_doEtaHighCut || _doEtaLowCut || _doBothEtaCut)
    cout << "   doEtaCut:  " << _theEtaLow << " < eta < " << _theEtaHigh << endl;
  if (_doAbsEtaHighCut || _doAbsEtaLowCut || _doBothAbsEtaCut)
    cout << "   doAbsEtaCut:  " << _theEtaLow << " < |eta| < " << _theEtaHigh << endl;
  if (_doPtHighCut || _doPtLowCut || _doBothPtCut)
    cout << "   doPtCut:  " << _thePtLow << " < pT < " << _thePtHigh << endl;
  if (_doPHighCut || _doPLowCut || _doBothPCut)
    cout << "   doPCut:  " << _thePLow << " < p < " << _thePHigh << endl;
  if (_doPzHighCut || _doPzLowCut || _doBothPzCut)
    cout << "   doPzCut:  " << _thePzLow << " < pz < " << _thePzHigh << endl;
  cout << "-----------------------------------------------------------------------" << endl;
}
