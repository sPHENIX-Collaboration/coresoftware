#include "ReadEICFiles.h"

#include "PHG4InEvent.h"
#include "PHG4Particlev1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/recoConsts.h>
#include <phool/getClass.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

// eicsmear classes
#include <eicsmear/erhic/EventMC.h>

// General Root and C++ classes
#include <TChain.h>

using namespace std;

///////////////////////////////////////////////////////////////////

ReadEICFiles::ReadEICFiles(const string &name):
  SubsysReco(name),
  filename(""),
  Tin(NULL),
  nEntries(0),
  entry(0),
  GenEvent(NULL)
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
  Tin = new TChain("EICTree", "EICTree");
  Tin->Add(name.c_str());
  GetTree();
  return true;
}

///////////////////////////////////////////////////////////////////

void ReadEICFiles::GetTree()
{
  /* Print the actual class of the event branch record,
     i.e. erhic::EventMilou or other */
  cout << "ReadEICFiles: Input Branch Event Class = "
       << Tin->GetBranch("event")->GetClassName() << endl;

  Tin->SetBranchAddress("event", &GenEvent);
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

  /* Check if there is an unused event left in input file */
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

  /* Get information about world geometry */
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

  /* Find InEvent node */
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve) {
    cout << PHWHERE << "no PHG4INEVENT node" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* Get event record from tree */
  Tin->GetEntry(entry);

  /* Loop over all particles for this event in input file */
  for (unsigned ii = 0; ii < GenEvent->GetNTracks(); ii++)
    {

      /* Get particle / track from event recors */
      erhic::VirtualParticle * track_ii = GenEvent->GetTrack(ii);

      /* Check if particle is stable final state particle */
      if ( track_ii->GetStatus()>10 )
        {
          continue;
        }

      /* get production vertex of particle */
      double vx = track_ii->GetVertex().X();
      double vy = track_ii->GetVertex().Y();
      double vz = track_ii->GetVertex().Z();
      double vt = 0;

      /* if world is a cylinder */
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
      /* if world is a box */
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
      /* if world is of unknown shape */
      else
        {
          cout << PHWHERE << " shape " << ishape << " not implemented. exiting" << endl;
          exit(1);
        }

      /* Add particle to Geant4 event */
      int idvtx = ineve->AddVtx(vx, vy, vz, vt);

      PHG4Particle *g4particle = new PHG4Particlev1();
      g4particle->set_pid(track_ii->Id());
      g4particle->set_px(track_ii->GetPx());
      g4particle->set_py(track_ii->GetPy());
      g4particle->set_pz(track_ii->GetPz());

      ineve->AddParticle(idvtx, g4particle);
    }

  /* Print event information if verbosity > 0 */
  if (verbosity > 0)
    {
      ineve->identify();
    }

  /* Count up number of 'used' events from input file */
  entry++;

  /* Done */
  return 0;
}
