#include "PHPy6ParticleTrigger.h"
#include "PHPy6GenTrigger.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>       // for GenVertex, GenVertex::particles_in...
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>     // for GenParticle
#include <HepMC/SimpleVector.h>    // for FourVector

#include <cmath>                  // for sqrt
#include <cstdlib>                // for abs
#include <iostream>

using namespace std;

//___________________________________________________________________________
PHPy6ParticleTrigger::PHPy6ParticleTrigger(const std::string &name)
  : PHPy6GenTrigger(name)
  ,

  _theEtaHigh(-999.9)
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

bool PHPy6ParticleTrigger::Apply(const HepMC::GenEvent *evt)
{
  // Print Out Trigger Information Once, for Posterity
  static int trig_info_printed = 0;
  if (trig_info_printed == 0)
  {
    PrintConfig();
    trig_info_printed = 1;
  }

  //  for ( HepMC::GenEvent::particle_const_iterator p
  //          = evt->particles_begin(); p != evt->particles_end(); ++p ){
  //    if ( (abs((*p)->pdg_id()) == 11) && ((*p)->status()==1) &&
  //         ((*p)->momentum().pseudoRapidity() > eta_low) && ((*p)->momentum().pseudoRapidity() < eta_high) &&
  //         ( ) {
  //      if(((*p)->pdg_id()) == 11) n_em_found++;
  //      if(((*p)->pdg_id()) == -11) n_ep_found++;
  //    }
  //  }
  //
  //  if( (RequireOR && ((n_em_found>=n_em_required)||(n_ep_found>=n_ep_required)) ) ||
  //      (RequireElectron && (n_em_found>=n_em_required)) ||
  //      (RequirePositron && (n_ep_found>=n_ep_required)) ||
  //      (RequireAND && (n_em_found>=n_em_required) && (n_ep_found>=n_ep_required)) ||
  //      (RequireCOMBO && (n_em_found+n_ep_found)>=n_comb_required) )
  //    {
  //      ++ntriggered_forward_electron;
  //      return true;
  //    }

  // Loop over all particles in the event
  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p)
  {
    // loop over all the trigger particle criteria
    for (int j = 0; j < int(_theParticles.size()); j++)
    {
      double p_pT = sqrt(pow((*p)->momentum().px(), 2) + pow((*p)->momentum().py(), 2));
      double p_pAbs = sqrt(pow((*p)->momentum().px(), 2) + pow((*p)->momentum().py(), 2) + pow((*p)->momentum().pz(), 2));
      if ((*p)->pdg_id() == _theParticles[j] &&
          (*p)->status() == 1)
      {  //only stable particles

        if (_doBothEtaCut && ((*p)->momentum().eta() < _theEtaLow ||
                              (*p)->momentum().eta() > _theEtaHigh)) continue;
        if (_doEtaLowCut && (*p)->momentum().eta() < _theEtaLow) continue;
        if (_doEtaHighCut && (*p)->momentum().eta() > _theEtaHigh) continue;

        if (_doBothAbsEtaCut && (abs((*p)->momentum().eta()) < _theEtaLow ||
                                 abs((*p)->momentum().eta()) > _theEtaHigh)) continue;
        if (_doAbsEtaLowCut && abs((*p)->momentum().eta()) < _theEtaLow) continue;
        if (_doAbsEtaHighCut && abs((*p)->momentum().eta()) > _theEtaHigh) continue;

        if (_doBothPtCut && (p_pT < _thePtLow ||
                             p_pT > _thePtHigh)) continue;
        if (_doPtHighCut && p_pT > _thePtHigh) continue;
        if (_doPtLowCut && p_pT < _thePtLow) continue;

        if (_doBothPCut && (p_pAbs < _thePLow ||
                            p_pAbs > _thePHigh)) continue;
        if (_doPHighCut && p_pAbs > _thePHigh) continue;
        if (_doPLowCut && p_pAbs < _thePLow) continue;

        if (_doBothPzCut && ((*p)->momentum().pz() < _thePzLow ||
                             (*p)->momentum().pz() > _thePzHigh)) continue;
        if (_doPzHighCut && (*p)->momentum().pz() > _thePzHigh) continue;
        if (_doPzLowCut && (*p)->momentum().pz() < _thePzLow) continue;

        if (Verbosity() > 5)
        {
          cout << "stable " << (*p)->pdg_id()
               << "  pt: " << p_pT
               << " pz: " << (*p)->momentum().pz()
               << " p: " << p_pAbs
               << " eta: " << (*p)->momentum().eta() << endl;
        }

        // loop over all partents to this particle
        bool passedParents = false;
        for (int k = 0; k < int(_theParents.size()); k++)
        {
          // check Mothers
          for (HepMC::GenVertex::particles_in_const_iterator p_parent = (*p)->production_vertex()->particles_in_const_begin();
               p_parent != (*p)->production_vertex()->particles_in_const_end();
               ++p_parent)
          {
            if (abs((*p_parent)->pdg_id()) == abs(_theParents[k]))
            {
              passedParents = true;
              if (Verbosity() > 5) cout << "found parent!" << endl;
              break;
            }
          }  //moms for loop
          if (passedParents) break;
        }  //parents for loop

        //If we made it here and it passes parents, success!
        if (_theParents.size() == 0 || passedParents) return true;

      }  //if _theParticles
    }    //_theParticles for loop

  }  //pythia event for loop

  return false;
}

void PHPy6ParticleTrigger::AddParticles(const std::string &particles)
{
  std::vector<int> addedParts = convertToInts(particles);
  _theParticles.insert(_theParticles.end(), addedParts.begin(), addedParts.end());
}

void PHPy6ParticleTrigger::AddParticles(int particle)
{
  _theParticles.push_back(particle);
}

void PHPy6ParticleTrigger::AddParticles(std::vector<int> particles)
{
  _theParticles.insert(_theParticles.end(), particles.begin(), particles.end());
}

void PHPy6ParticleTrigger::AddParents(const std::string &parents)
{
  std::vector<int> addedParents = convertToInts(parents);
  _theParents.insert(_theParents.end(), addedParents.begin(), addedParents.end());
}

void PHPy6ParticleTrigger::AddParents(int parent)
{
  _theParents.push_back(parent);
}

void PHPy6ParticleTrigger::AddParents(std::vector<int> parents)
{
  _theParents.insert(_theParents.end(), parents.begin(), parents.end());
}

void PHPy6ParticleTrigger::SetPtHigh(double pt)
{
  _thePtHigh = pt;
  if (_doPtLowCut)
    _doBothPtCut = true;
  else
    _doPtHighCut = true;
}

void PHPy6ParticleTrigger::SetPtLow(double pt)
{
  _thePtLow = pt;
  if (_doPtHighCut)
    _doBothPtCut = true;
  else
    _doPtLowCut = true;
}

void PHPy6ParticleTrigger::SetPtHighLow(double ptHigh, double ptLow)
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

void PHPy6ParticleTrigger::SetPHigh(double p)
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

void PHPy6ParticleTrigger::SetPLow(double p)
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

void PHPy6ParticleTrigger::SetPHighLow(double pHigh, double pLow)
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

void PHPy6ParticleTrigger::SetEtaHigh(double eta)
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

void PHPy6ParticleTrigger::SetEtaLow(double eta)
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

void PHPy6ParticleTrigger::SetEtaHighLow(double etaHigh, double etaLow)
{
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;
  _doBothEtaCut = true;
  _doEtaHighCut = false;
  _doEtaLowCut = false;
}

void PHPy6ParticleTrigger::SetAbsEtaHigh(double eta)
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

void PHPy6ParticleTrigger::SetAbsEtaLow(double eta)
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

void PHPy6ParticleTrigger::SetAbsEtaHighLow(double etaHigh, double etaLow)
{
  _theEtaHigh = etaHigh;
  _theEtaLow = etaLow;
  _doBothAbsEtaCut = true;
  _doAbsEtaLowCut = false;
  _doAbsEtaHighCut = false;
}

void PHPy6ParticleTrigger::SetPzHigh(double pz)
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

void PHPy6ParticleTrigger::SetPzLow(double pz)
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

void PHPy6ParticleTrigger::SetPzHighLow(double pzHigh, double pzLow)
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

void PHPy6ParticleTrigger::PrintConfig()
{
  cout << "---------------- PHPy6ParticleTrigger::PrintConfig --------------------" << endl;
  cout << "   Particles: ";
  for (int i = 0; i < int(_theParticles.size()); i++) cout << _theParticles[i] << "  ";
  cout << endl;

  cout << "   Parents: ";
  for (int i = 0; i < int(_theParents.size()); i++) cout << _theParents[i] << "  ";
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
