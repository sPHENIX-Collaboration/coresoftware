#include <PHG4PhPyReader.h>
#include <PHPythiaContainer.h>
#include <PHPyCommon.h>

#include <getClass.h>
#include <Fun4AllReturnCodes.h>

#include <TLorentzVector.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8) 
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include "PHG4InEvent.h"
#include "PHG4Particlev1.h"

#include <Fun4AllReturnCodes.h>
#include <recoConsts.h>

#include <PHCompositeNode.h>
#include <PHNodeIterator.h>

#include <getClass.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <CLHEP/Random/RandGauss.h>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4PhysicalConstants.hh>

using namespace std;

PHG4PhPyReader::PHG4PhPyReader(const std::string &name) :
    SubsysReco(name), 
    _input_node("PHPythia"),
    phpythia(NULL), 
    _event_vx(0), 
    _event_vy(0), 
    _event_vz(0)
{
  _vtx_offset.resize(3, 0);
  _vtx_sigma.resize(3, 0);
}

PHG4PhPyReader::~PHG4PhPyReader()
{
}

int
PHG4PhPyReader::Init(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

int
PHG4PhPyReader::End(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

int
PHG4PhPyReader::process_event(PHCompositeNode *topNode)
{

  recoConsts *rc = recoConsts::instance();
  float worldsizex = rc->get_FloatFlag("WorldSizex");
  float worldsizey = rc->get_FloatFlag("WorldSizey");
  float worldsizez = rc->get_FloatFlag("WorldSizez");
  string worldshape = rc->get_CharFlag("WorldShape");
  enum
  {
    ShapeG4Tubs = 0, ShapeG4Box = 1
  };
  int ishape;
  if (worldshape == "G4Tubs")
    {
      ishape = ShapeG4Tubs;
    }
  else if (worldshape == "G4Box")
    {
      ishape = ShapeG4Box;
    }
  else
    {
      cout << PHWHERE << " unknown world shape " << worldshape << endl;
      exit(1);
    }
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");

  // Get PYTHIA Particles
  phpythia = findNode::getClass<PHPythiaContainer>(topNode, _input_node.c_str());
  if (!phpythia)
    {
      cout << PHWHERE << "Unable to get " << _input_node << ", is Node missing?" << endl;
      return ABORTEVENT;
    }

  int npart = phpythia->size();
  if (verbosity)
    cout << "PHG4PhPyReader::process_event - in put list size: " << npart
        << endl;

  // vertex

  _event_vx = CLHEP::RandGauss::shoot(_vtx_offset[0], _vtx_sigma[0]);
  _event_vy = CLHEP::RandGauss::shoot(_vtx_offset[1], _vtx_sigma[1]);
  _event_vz = CLHEP::RandGauss::shoot(_vtx_offset[2], _vtx_sigma[2]);

  if (verbosity)
    cout << "PHG4PhPyReader::process_event - event vertex = " << _event_vx
        << ", " << _event_vy << ", " << _event_vz << "" << endl;

  for (int ipart = 0; ipart < npart; ipart++)
    {
      TMCParticle *Particle = phpythia->getParticle(ipart);

      //  Unit convension for TMCParticle
      //      Float_t fEnergy Energy [GeV] ( LUJETS P[4] )
      //      Int_t fFirstChild id of first child ( LUJETS K[4] )
      //      Int_t fKF KF flavour code ( LUJETS K[2] )
      //      Int_t fKS status of particle ( LUJETS K[1] )
      //      Int_t fLastChild  id of last child ( LUJETS K[5] )
      //      Float_t fLifetime proper lifetime [mm/c] ( LUJETS V[5] )
      //      Float_t fMass Mass [Gev/c^2] ( LUJETS P[5] )
      //      Int_t fParent parrent's id ( LUJETS K[3] )
      //      Float_t fPx X momenta [GeV/c] ( LUJETS P[1] )
      //      Float_t fPy Y momenta [GeV/c] ( LUJETS P[2] )
      //      Float_t fPz Z momenta [GeV/c] ( LUJETS P[3] )
      //      Float_t fTime time of procuction [mm/c]( LUJETS V[4] )
      //      Float_t fVx X vertex [mm] ( LUJETS V[1] )
      //      Float_t fVy Y vertex [mm] ( LUJETS V[2] )
      //      Float_t fVz Z vertex [mm] ( LUJETS V[3] )


//      if (Particle->GetFirstChild())
      if (Particle->GetKS() < kMIN_ALIVE_KS
          or Particle->GetKS() > kMAX_ALIVE_KS)
        {
          if (Verbosity() >= 2)
            {
              cout << "PHG4PhPyReader::process_event - drop decayed particle "
                  << Particle->GetName() << " status " << Particle->GetKS()
                  << endl;
            }

          continue;
        }

      // get production vertex of particle
      double vx = Particle->GetVx()*mm/cm  + _event_vx;
      double vy = Particle->GetVy()*mm/cm + _event_vy;
      double vz = Particle->GetVz()*mm/cm + _event_vz;
      double vt = Particle->GetTime()*(mm/c_light)/s;
      // world is a cylinder
      if (ishape == ShapeG4Tubs)
        {
          if (sqrt(vx * vx + vy * vy) > worldsizey / 2
              || fabs(vz) > worldsizez / 2)
            {
              cout << "vertex x/y/z" << vx << "/" << vy << "/" << vz
                  << " outside world volume radius/z (+-) " << worldsizex / 2
                  << "/" << worldsizez / 2 << ", dropping it and its particles"
                  << endl;
              continue;
            }
        }
      else if (ishape == ShapeG4Box)
        {
          if (fabs(vx) > worldsizex / 2 || fabs(vy) > worldsizey / 2
              || fabs(vz) > worldsizez / 2)
            {
              cout << "vertex x/y/z" << vx << "/" << vy << "/" << vz
                  << " outside world volume x/y/z (+-) " << worldsizex / 2
                  << "/" << worldsizey / 2 << "/" << worldsizez / 2
                  << ", dropping it and its particles" << endl;
              continue;
            }
        }
      else
        {
          cout << PHWHERE << " shape " << ishape << " not implemented. exiting"
              << endl;
          exit(1);
        }
      int idvtx = ineve->AddVtx(vx, vy, vz, vt);
      PHG4Particle *g4particle = new PHG4Particlev1();
      g4particle->set_pid(Particle->GetKF());
      g4particle->set_px(Particle->GetPx());
      g4particle->set_py(Particle->GetPy());
      g4particle->set_pz(Particle->GetPz());
      g4particle->set_name(Particle->GetName());
      ineve->AddParticle(idvtx, g4particle);

      if (Verbosity() >= 2)
        {

          cout << "PHG4PhPyReader::process_event - generate particle "
              << Particle->GetName() << " @ vertex = " << vx << ", " << vy
              << ", " << vz << " cm, t = " << vt << " s, momentum = "
              << Particle->GetPx() << ", " << Particle->GetPy() << ", "
              << Particle->GetPz() <<" GeV/c" << endl;
        }
    }
  if (verbosity > 0)
    {
      ineve->identify();
    }

  return EVENT_OK; // trigger condition not found, don't write out event
}

