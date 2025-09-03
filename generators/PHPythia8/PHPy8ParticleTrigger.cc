#include "PHPy8ParticleTrigger.h"

#include <Pythia8/Event.h>  // for Event, Particle
#include <Pythia8/Pythia.h>

#include <algorithm>  // for max
#include <cstdlib>    // for abs
#include <iostream>   // for operator<<, endl, basic_ostream, basic_o...

//__________________________________________________________
PHPy8ParticleTrigger::PHPy8ParticleTrigger(const std::string &name)
  : PHPy8GenTrigger(name)
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

PHPy8ParticleTrigger::~PHPy8ParticleTrigger()
{
  if (Verbosity() > 0)
  {
    PrintConfig();
  }
}

bool PHPy8ParticleTrigger::Apply(Pythia8::Pythia *pythia)
{
  if (Verbosity() > 2)
  {
    std::cout << "PHPy8ParticleTrigger::Apply - pythia event size: "
         << pythia->event.size() << std::endl;
  }

  // Loop over all particles in the event
  for (int i = 0; i < pythia->event.size(); ++i)
  {
    // loop over all the trigger particle criteria
    for (int _theParticle : _theParticles)
    {
      if (pythia->event[i].id() == _theParticle &&
          (pythia->event[i].status() > 0    // only stable particles
           || (! m_doStableParticleOnly)  // or not
           ))
      {
        if (_doBothYCut && (pythia->event[i].y() < _theYLow ||
                            pythia->event[i].y() > _theYHigh))
        {
          continue;
        }
        if (_doYLowCut && pythia->event[i].y() < _theYLow)
        {
          continue;
        }
        if (_doYHighCut && pythia->event[i].y() > _theYHigh)
        {
          continue;
        }

        if (_doBothEtaCut && (pythia->event[i].eta() < _theEtaLow ||
                              pythia->event[i].eta() > _theEtaHigh))
        {
          continue;
        }
        if (_doEtaLowCut && pythia->event[i].eta() < _theEtaLow)
        {
          continue;
        }
        if (_doEtaHighCut && pythia->event[i].eta() > _theEtaHigh)
        {
          continue;
        }

        if (_doBothAbsEtaCut && (std::abs(pythia->event[i].eta()) < _theEtaLow ||
                                 std::abs(pythia->event[i].eta()) > _theEtaHigh))
        {
          continue;
        }
        if (_doAbsEtaLowCut && std::abs(pythia->event[i].eta()) < _theEtaLow)
        {
          continue;
        }
        if (_doAbsEtaHighCut && std::abs(pythia->event[i].eta()) > _theEtaHigh)
        {
          continue;
        }

        if (_doBothPtCut && (pythia->event[i].pT() < _thePtLow ||
                             pythia->event[i].pT() > _thePtHigh))
        {
          continue;
        }
        if (_doPtHighCut && pythia->event[i].pT() > _thePtHigh)
        {
          continue;
        }
        if (_doPtLowCut && pythia->event[i].pT() < _thePtLow)
        {
          continue;
        }

        if (_doBothPCut && (pythia->event[i].pAbs() < _thePLow ||
                            pythia->event[i].pAbs() > _thePHigh))
        {
          continue;
        }
        if (_doPHighCut && pythia->event[i].pAbs() > _thePHigh)
        {
          continue;
        }
        if (_doPLowCut && pythia->event[i].pAbs() < _thePLow)
        {
          continue;
        }

        if (_doBothPzCut && (pythia->event[i].pz() < _thePzLow ||
                             pythia->event[i].pz() > _thePzHigh))
        {
          continue;
        }
        if (_doPzHighCut && pythia->event[i].pz() > _thePzHigh)
        {
          continue;
        }
        if (_doPzLowCut && pythia->event[i].pz() < _thePzLow)
        {
          continue;
        }

        if (Verbosity() > 5)
        {
          std::cout << "stable " << pythia->event[i].id()
               << "  pt: " << pythia->event[i].pT()
               << " pz: " << pythia->event[i].pz()
               << " p: " << pythia->event[i].pAbs()
               << " eta: " << pythia->event[i].eta()
               << " y: " << pythia->event[i].y() << std::endl;
        }

        // loop over all partents to this particle
        bool passedParents = false;
        for (int _theParent : _theParents)
        {
          // check Mothers
          std::vector<int> moms = pythia->event[i].motherList();
          for (int mom : moms)
          {
            if (std::abs(pythia->event[mom].id()) == std::abs(_theParent))
            {
              passedParents = true;
              if (Verbosity() > 5)
              {
                std::cout << "found parent!" << std::endl;
              }
              break;
            }
          }  // moms for loop
          if (passedParents)
          {
            break;
          }
        }  // parents for loop

        // If we made it here and it passes parents, success!
        if (_theParents.size() == 0 || passedParents)
        {
          return true;
        }

      }  // if _theParticles
    }    //_theParticles for loop

  }  // pythia event for loop

  return false;
}

void PHPy8ParticleTrigger::AddParticles(const std::string &particles)
{
  std::vector<int> addedParts = convertToInts(particles);
  _theParticles.insert(_theParticles.end(), addedParts.begin(), addedParts.end());
}

void PHPy8ParticleTrigger::AddParticles(int particle)
{
  _theParticles.push_back(particle);
}

void PHPy8ParticleTrigger::AddParticles(std::vector<int> particles)
{
  _theParticles.insert(_theParticles.end(), particles.begin(), particles.end());
}

void PHPy8ParticleTrigger::AddParents(const std::string &parents)
{
  std::vector<int> addedParents = convertToInts(parents);
  _theParents.insert(_theParents.end(), addedParents.begin(), addedParents.end());
}

void PHPy8ParticleTrigger::AddParents(int parent)
{
  _theParents.push_back(parent);
}

void PHPy8ParticleTrigger::AddParents(std::vector<int> parents)
{
  _theParents.insert(_theParents.end(), parents.begin(), parents.end());
}

void PHPy8ParticleTrigger::SetPtHigh(double pt)
{
  _thePtHigh = pt;
  if (_doPtLowCut)
  {
    _doBothPtCut = true;
  }
  else
  {
    _doPtHighCut = true;
  }
}

void PHPy8ParticleTrigger::SetPtLow(double pt)
{
  _thePtLow = pt;
  if (_doPtHighCut)
  {
    _doBothPtCut = true;
  }
  else
  {
    _doPtLowCut = true;
  }
}

void PHPy8ParticleTrigger::SetPtHighLow(double ptHigh, double ptLow)
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

void PHPy8ParticleTrigger::SetPHigh(double p)
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

void PHPy8ParticleTrigger::SetPLow(double p)
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

void PHPy8ParticleTrigger::SetPHighLow(double pHigh, double pLow)
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

void PHPy8ParticleTrigger::SetYHigh(double Y)
{
  _theYHigh = Y;
  if (_doYLowCut)
  {
    _doBothYCut = true;
    _doYLowCut = false;
  }
  else
  {
    _doYHighCut = true;
  }
}

void PHPy8ParticleTrigger::SetYLow(double Y)
{
  _theYLow = Y;
  if (_doYHighCut)
  {
    _doBothYCut = true;
    _doYHighCut = false;
  }
  else
  {
    _doYLowCut = true;
  }
}

void PHPy8ParticleTrigger::SetYHighLow(double YHigh, double YLow)
{
  _theYHigh = YHigh;
  _theYLow = YLow;
  _doBothYCut = true;
  _doYHighCut = false;
  _doYLowCut = false;
}

void PHPy8ParticleTrigger::SetEtaHigh(double eta)
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

void PHPy8ParticleTrigger::SetEtaLow(double eta)
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

void PHPy8ParticleTrigger::SetEtaHighLow(double etaHigh, double etaLow)
{
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;
  _doBothEtaCut = true;
  _doEtaHighCut = false;
  _doEtaLowCut = false;
}

void PHPy8ParticleTrigger::SetAbsEtaHigh(double eta)
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

void PHPy8ParticleTrigger::SetAbsEtaLow(double eta)
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

void PHPy8ParticleTrigger::SetAbsEtaHighLow(double etaHigh, double etaLow)
{
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;
  _doBothAbsEtaCut = true;
  _doAbsEtaLowCut = false;
  _doAbsEtaHighCut = false;
}

void PHPy8ParticleTrigger::SetPzHigh(double pz)
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

void PHPy8ParticleTrigger::SetPzLow(double pz)
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

void PHPy8ParticleTrigger::SetPzHighLow(double pzHigh, double pzLow)
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

void PHPy8ParticleTrigger::PrintConfig()
{
  std::cout << "---------------- PHPy8ParticleTrigger::PrintConfig --------------------" << std::endl;

  if (m_doStableParticleOnly)
  {
    std::cout << "Process stable particles only." << std::endl;
  }
  else
  {
    std::cout << "Process both unstable and stable particles." << std::endl;
  }

  std::cout << "   Particles: ";
  for (int _theParticle : _theParticles)
  {
    std::cout << _theParticle << "  ";
  }
  std::cout << std::endl;

  std::cout << "   Parents: ";
  for (int _theParent : _theParents)
  {
    std::cout << _theParent << "  ";
  }
  std::cout << std::endl;

  if (_doYHighCut || _doYLowCut || _doBothYCut)
  {
    std::cout << "   doYCut:  " << _theYLow << " < Y < " << _theYHigh << std::endl;
  }
  if (_doEtaHighCut || _doEtaLowCut || _doBothEtaCut)
  {
    std::cout << "   doEtaCut:  " << _theEtaLow << " < eta < " << _theEtaHigh << std::endl;
  }
  if (_doAbsEtaHighCut || _doAbsEtaLowCut || _doBothAbsEtaCut)
  {
    std::cout << "   doAbsEtaCut:  " << _theEtaLow << " < |eta| < " << _theEtaHigh << std::endl;
  }
  if (_doPtHighCut || _doPtLowCut || _doBothPtCut)
  {
    std::cout << "   doPtCut:  " << _thePtLow << " < pT < " << _thePtHigh << std::endl;
  }
  if (_doPHighCut || _doPLowCut || _doBothPCut)
  {
    std::cout << "   doPCut:  " << _thePLow << " < p < " << _thePHigh << std::endl;
  }
  if (_doPzHighCut || _doPzLowCut || _doBothPzCut)
  {
    std::cout << "   doPzCut:  " << _thePzLow << " < pz < " << _thePzHigh << std::endl;
  }
  std::cout << "-----------------------------------------------------------------------" << std::endl;
}
