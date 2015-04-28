#include "ReadEICFiles.h"

#include "PHG4InEvent.h"
#include "PHG4Particlev1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/recoConsts.h>
#include <fun4all/getClass.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>


//  General Root, Pythia and C classes
#include <TParticle.h>
#include <TClonesArray.h>

#include <TChain.h>

#include <vector>
#include <fstream>

using namespace std;

///////////////////////////////////////////////////////////////////

ReadEICFiles::ReadEICFiles(const string &name):
  SubsysReco(name),
  Tin(NULL),
  Particle(NULL),
  ParticleArray(NULL),
  nEntries(0),
  entry(0)
{
  return;
}

///////////////////////////////////////////////////////////////////

ReadEICFiles::~ReadEICFiles()
{
  delete Tin;
}

///////////////////////////////////////////////////////////////////


bool
ReadEICFiles::OpenInputFile(const string &name)
{
  filename = name;
  Tin = new TChain("T", "T");
  Tin->Add(name.c_str());
  GetTree();
  return true;
}

///////////////////////////////////////////////////////////////////

void ReadEICFiles::GetTree()
{
  Tin->SetBranchAddress("process_id", &ProcessID);
  Tin->SetBranchAddress("y", &Y);
  Tin->SetBranchAddress("Q2", &Q2);
  Tin->SetBranchAddress("x", &X);
  Tin->SetBranchAddress("nu", &NU);
  Tin->SetBranchAddress("W2", &W2);
  Tin->SetBranchAddress("particles", &ParticleArray);
  nEntries = Tin->GetEntries();
}

int
ReadEICFiles::Init(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
    {
      PHNodeIterator iter( topNode );
      PHCompositeNode *dstNode;
      dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

      ineve = new PHG4InEvent();
      PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
      dstNode->addNode(newNode);
    }
  return 0;
}


int
ReadEICFiles::process_event(PHCompositeNode *topNode)
{
  if (entry >= nEntries)
    {
      if (filename.size() > 0)
        {
          cout << "input file " << filename << " exhausted" << endl;
        }
      else
        {
          cout << "no file opened" << endl;
        }
      return Fun4AllReturnCodes::ABORTRUN;
    }

  recoConsts *rc = recoConsts::instance();
  float worldsizex = rc->get_FloatFlag("WorldSizex");
  float worldsizey = rc->get_FloatFlag("WorldSizey");
  float worldsizez = rc->get_FloatFlag("WorldSizez");
  string worldshape = rc->get_CharFlag("WorldShape");
  enum {ShapeG4Tubs = 0, ShapeG4Box = 1};
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

  Tin->GetEntry(entry);
  for (int ii = 0; ii < ParticleArray->GetEntries(); ii++)
    {
      Particle = (TParticle *)ParticleArray->At(ii);
      //skip particles which have daughters
      //      if (Particle->GetFirstDaughter())
      if (Particle->GetStatusCode()>10)
        {
          continue;
        }
      // get production vertex of particle
      double vx = Particle->Vx();
      double vy = Particle->Vy();
      double vz = Particle->Vz();
      double vt = Particle->T();
      // world is a cylinder
      if (ishape == ShapeG4Tubs)
        {
          if (sqrt(vx*vx + vy*vy) > worldsizey / 2 || fabs(vz) > worldsizez / 2)
            {
              cout << "vertex x/y/z" << vx << "/" << vy << "/" << vz
                   << " outside world volume radius/z (+-) " << worldsizex / 2 << "/"
                   << worldsizez / 2
                   << ", dropping it and its particles" << endl;
              continue;
            }
        }
      else if (ishape == ShapeG4Box)
        {
          if (fabs(vx) > worldsizex / 2 || fabs(vy) > worldsizey / 2 || fabs(vz) > worldsizez / 2)
            {
              cout << "vertex x/y/z" << vx << "/" << vy << "/" << vz
                   << " outside world volume x/y/z (+-) " << worldsizex / 2 << "/"
                   << worldsizey / 2 << "/" << worldsizez / 2
                   << ", dropping it and its particles" << endl;
              continue;
            }
        }
      else
        {
          cout << PHWHERE << " shape " << ishape << " not implemented. exiting" << endl;
          exit(1);
        }
      int idvtx = ineve->AddVtx(vx, vy, vz, vt);
      PHG4Particle *g4particle = new PHG4Particlev1();
      g4particle->set_pid(Particle->GetPdgCode());
      g4particle->set_px(Particle->Px());
      g4particle->set_py(Particle->Py());
      g4particle->set_pz(Particle->Pz());
      ineve->AddParticle(idvtx, g4particle);
    }
  if (verbosity > 0)
    {
      ineve->identify();
    }
  entry++;
  return 0;
}
