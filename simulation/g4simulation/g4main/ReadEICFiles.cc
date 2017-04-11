#include "ReadEICFiles.h"

#include <phhepmc/PHHepMCGenEvent.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>

// eicsmear classes
#include <eicsmear/erhic/EventMC.h>

// General Root and C++ classes
#include <TChain.h>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;

///////////////////////////////////////////////////////////////////

ReadEICFiles::ReadEICFiles(const string &name):
  SubsysReco(name),
  filename(""),
  Tin(NULL),
  nEntries(0),
  entry(0),
  GenEvent(NULL),
  _node_name("PHHepMCGenEvent")
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
  /* Create node tree */
  CreateNodeTree(topNode);
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

  /* Get event record from input file */
  Tin->GetEntry(entry);

  /* Create GenEvent */
  HepMC::GenEvent* evt = new HepMC::GenEvent();

  /* define the units (Pythia uses GeV and mm) */
  evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

  /* add global information to the event */
  evt->set_event_number(entry);

  // -----------------



  /* Create dummy vertex at 0 for HepMC event */
  HepMC::GenVertex* hepmcvtx = new HepMC::GenVertex(HepMC::FourVector(0,0,0,0));

  /* Loop over all particles for this event in input file */
  for (unsigned ii = 0; ii < GenEvent->GetNTracks(); ii++)
    {
      /* Get particle / track from event records.
       * Class documentation for erhic::VirtualParticle at
       * http://www.star.bnl.gov/~tpb/eic-smear/classerhic_1_1_virtual_particle.html */
      erhic::VirtualParticle * track_ii = GenEvent->GetTrack(ii);

      cout << "index: " << ii <<" , ID: " << track_ii->Id() << " , parent index: " << track_ii->GetParentIndex() << " , vertex: (" << track_ii->GetVertex()[0] << ", " <<  track_ii->GetVertex()[1] << ", " <<  track_ii->GetVertex()[2] << ") , status:" << track_ii->GetStatus() << endl;

      /* Create HepMC particle record */
      HepMC::GenParticle *hepmcpart = new HepMC::GenParticle( HepMC::FourVector(track_ii->GetPx(),
										track_ii->GetPy(),
										track_ii->GetPz(),
										track_ii->GetE()),
							      track_ii->Id());

      /* translate eic-smear status codes to HepMC status codes */
      switch ( track_ii->GetStatus() )
	{
	case 1: hepmcpart->set_status( 1 ); break;  // 'stable particle'

	case 21: hepmcpart->set_status( 3 ); break; // 'documentation line'

	default: hepmcpart->set_status( 0 ); break; // 'null entry'
	}

      /* if ParentIndex == 0, then this is a beam particle (which gets HepMC status 4)*/
      if ( track_ii->GetParentIndex() == 0 )
	hepmcpart->set_status( 4 );

      /* add particle to Hep MC vertex */
      hepmcvtx->add_particle_out(hepmcpart);
    }

  /* Add HepMC dummy vertex to event */
  evt->add_vertex(hepmcvtx);

  /* pass HepMC to PHNode*/
  bool success = _phhepmcevt->addEvent(evt);
  if (!success) {
    cout << "ReadEICFiles::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /* Count up number of 'used' events from input file */
  entry++;

  /* Done */
  return 0;
}

int ReadEICFiles::CreateNodeTree(PHCompositeNode *topNode) {

  /* HepMC node */
  PHCompositeNode *dstNode;
  PHNodeIterator iter(topNode);

  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing doing nothing" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _phhepmcevt = new PHHepMCGenEvent();
  PHObjectNode_t *newNode = new PHObjectNode_t(_phhepmcevt,_node_name.c_str(),"PHObject");
  dstNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}
