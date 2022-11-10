#include "ReadEICFiles.h"

#include "EicEventHeader.h"
#include "EicEventHeaderv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>  // for GenParticle
#include <HepMC/PdfInfo.h>       // for PdfInfo
#include <HepMC/SimpleVector.h>  // for FourVector
#include <HepMC/Units.h>         // for GEV, MM

// eicsmear classes
#include <eicsmear/erhic/EventMC.h>
#include <eicsmear/erhic/EventMilou.h>
#include <eicsmear/erhic/ParticleMC.h>  // for ParticleMC
#include <eicsmear/erhic/Pid.h>         // for Pid

#include <eicsmear/erhic/EventDEMP.h>

// General Root and C++ classes
#include <TBranch.h>  // for TBranch
#include <TChain.h>
#include <TVector3.h>  // for TVector3

#include <iostream>  // for operator<<, endl, basic_ostream
#include <vector>    // for vector

class PHHepMCGenEvent;

using namespace std;

///////////////////////////////////////////////////////////////////

ReadEICFiles::ReadEICFiles(const string &name)
  : SubsysReco(name)
  , filename("")
  , Tin(nullptr)
  , nEntries(0)
  , entry(0)
  , m_EvtGenId(EvtGen::Unknown)
  , GenEvent(nullptr)
  , _node_name("PHHepMCGenEvent")
{
  PHHepMCGenHelper::set_embedding_id(1);  // default embedding ID to 1
  return;
}

///////////////////////////////////////////////////////////////////

ReadEICFiles::~ReadEICFiles()
{
  delete Tin;
}

///////////////////////////////////////////////////////////////////

bool ReadEICFiles::OpenInputFile(const string &name)
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
  string EventClass(Tin->GetBranch("event")->GetClassName());
  if (Verbosity() > 0)
  {
    cout << "ReadEICFiles: Input Branch Event Class = "
         << EventClass << endl;
  }
  if (EventClass.find("Milou") != string::npos)
  {
    m_EvtGenId = EvtGen::Milou;
  }

  if (EventClass.find("DEMP") != string::npos)
  {
    m_EvtGenId = EvtGen::DEMP;
  }

  Tin->SetBranchAddress("event", &GenEvent);
  nEntries = Tin->GetEntries();
}

int ReadEICFiles::Init(PHCompositeNode *topNode)
{
  /* Create node tree */
  CreateNodeTree(topNode);
  return 0;
}

int ReadEICFiles::process_event(PHCompositeNode *topNode)
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

  EicEventHeader *evthead = findNode::getClass<EicEventHeader>(topNode, "EicEventHeader");
  switch (m_EvtGenId)
  {
  case EvtGen::Milou:
  {
    erhic::EventMilou *gen = dynamic_cast<erhic::EventMilou *>(GenEvent);
    evthead->set_eventgenerator_type(EicEventHeader::EvtGen::Milou);
    evthead->set_milou_weight(gen->weight);
    evthead->set_milou_trueX(gen->trueX);
    evthead->set_milou_trueQ2(gen->trueQ2);
  }
  break;
  case EvtGen::DEMP:
  {
    erhic::EventDEMP *gen = dynamic_cast<erhic::EventDEMP *>(GenEvent);
    evthead->set_eventgenerator_type(EicEventHeader::EvtGen::DEMP);
    evthead->set_demp_weight(gen->weight);
  }
  break;
  case EvtGen::Unknown:
    cout << "unknown event generator" << endl;
    break;
  default:
    cout << "what is this " << m_EvtGenId << " ????" << endl;
    break;
  }

  /* Create GenEvent */
  HepMC::GenEvent *evt = new HepMC::GenEvent();

  /* define the units (Pythia uses GeV and mm) */
  evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

  /* add global information to the event */
  evt->set_event_number(entry);

  /* process ID from pythia */
  evt->set_signal_process_id(GenEvent->GetProcess());

  /* Set the PDF information */
  HepMC::PdfInfo pdfinfo;
  pdfinfo.set_x1(1);
  pdfinfo.set_x2(GenEvent->GetX());
  pdfinfo.set_scalePDF(GenEvent->GetQ2());
  evt->set_pdf_info(pdfinfo);

  /* Loop over all particles for this event in input file and fill
   * vector with HepMC particles */
  vector<HepMC::GenParticle *> hepmc_particles;
  vector<unsigned> origin_index;

  /* save pointers to beam particles */
  HepMC::GenParticle *hepmc_beam1 = nullptr;
  HepMC::GenParticle *hepmc_beam2 = nullptr;

  for (unsigned ii = 0; ii < GenEvent->GetNTracks(); ii++)
  {
    /* Get particle / track from event records.
     * Class documentation for erhic::VirtualParticle at
     * http://www.star.bnl.gov/~tpb/eic-smear/classerhic_1_1_virtual_particle.html */
    erhic::ParticleMC *track_ii = GenEvent->GetTrack(ii);

    if (Verbosity() > 1)
    {
      cout << __PRETTY_FUNCTION__ << " : " << __LINE__ << endl;
      cout << "\ttrack " << ii << endl;
      cout << "\t4mom\t" << track_ii->GetPx() << "\t" << track_ii->GetPy() << "\t" << track_ii->GetPz() << "\t" << track_ii->GetE() << endl;
      cout << "\tstatus= " << track_ii->GetStatus() << "\tindex= " << track_ii->GetIndex() << "\tmass= " << track_ii->GetM() << endl;
    }

    /* Create HepMC particle record */
    HepMC::GenParticle *hepmcpart = new HepMC::GenParticle(HepMC::FourVector(track_ii->GetPx(),
                                                                             track_ii->GetPy(),
                                                                             track_ii->GetPz(),
                                                                             track_ii->GetE()),
                                                           track_ii->Id());

    /* translate eic-smear status codes to HepMC status codes */
    switch (track_ii->GetStatus())
    {
    case 1:
      hepmcpart->set_status(1);
      break;  // 'stable particle'

    case 21:
      hepmcpart->set_status(3);
      break;  // 'documentation line'

    default:
      hepmcpart->set_status(0);
      break;  // 'null entry'
    }

    /* assume the first two particles are the beam particles (which getsHepMC status 4)*/
    if (ii < 2)
      hepmcpart->set_status(4);

    /* add particle information */
    hepmcpart->setGeneratedMass(track_ii->GetM());

    /* append particle to vector */
    hepmc_particles.push_back(hepmcpart);
    origin_index.push_back(track_ii->GetIndex());

    /* if first particle, call this the first beam particle */
    if (ii == 0)
      hepmc_beam1 = hepmcpart;

    /* if second particle, call this the second beam particle */
    if (ii == 1)
      hepmc_beam2 = hepmcpart;
  }

  /* Check if hepmc_particles and origin_index vectors are the same size */
  if (hepmc_particles.size() != origin_index.size())
  {
    cout << "ReadEICFiles::process_event - Lengths of HepMC particles and Origin index vectors do not match!" << endl;

    delete evt;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /* add HepMC particles to Hep MC vertices; skip first two particles
   * in loop, assuming that they are the beam particles */
  vector<HepMC::GenVertex *> hepmc_vertices;

  for (unsigned p = 2; p < hepmc_particles.size(); p++)
  {
    HepMC::GenParticle *pp = hepmc_particles.at(p);

    /* continue if vertices for particle are already set */
    if (pp->production_vertex() && pp->end_vertex())
      continue;

    /* access mother particle vertex */
    erhic::ParticleMC *track_pp = GenEvent->GetTrack(p);

    unsigned parent_index = track_pp->GetParentIndex();

    HepMC::GenParticle *pmother = nullptr;
    for (unsigned m = 0; m < hepmc_particles.size(); m++)
    {
      if (origin_index.at(m) == parent_index)
      {
        pmother = hepmc_particles.at(m);
        break;
      }
    }

    /* if mother does not exist: create new vertex and add this particle as outgoing to vertex */
    if (!pmother)
    {
      HepMC::GenVertex *hepmcvtx = new HepMC::GenVertex(HepMC::FourVector(track_pp->GetVertex()[0],
                                                                          track_pp->GetVertex()[1],
                                                                          track_pp->GetVertex()[2],
                                                                          0));
      hepmc_vertices.push_back(hepmcvtx);
      hepmcvtx->add_particle_out(pp);
      continue;
    }
    /* if mother exists and has end vertex: add this particle as outgoing to the mother's end vertex */
    else if (pmother->end_vertex())
    {
      pmother->end_vertex()->add_particle_out(pp);
    }
    /* if mother exists and has no end vertex: create new vertex */
    else
    {
      HepMC::GenVertex *hepmcvtx = new HepMC::GenVertex(HepMC::FourVector(track_pp->GetVertex()[0],
                                                                          track_pp->GetVertex()[1],
                                                                          track_pp->GetVertex()[2],
                                                                          0));
      hepmc_vertices.push_back(hepmcvtx);
      hepmcvtx->add_particle_in(pmother);
      pmother->end_vertex()->add_particle_out(pp);
    }
  }

  /* Add end vertex to beam particles if they don't have one yet */
  for (unsigned p = 0; p < 2; p++)
  {
    HepMC::GenParticle *pp = hepmc_particles.at(p);

    if (!pp->end_vertex())
    {
      /* create collision vertex */
      HepMC::GenVertex *hepmcvtx = new HepMC::GenVertex(HepMC::FourVector(0,
                                                                          0,
                                                                          0,
                                                                          0));
      hepmc_vertices.push_back(hepmcvtx);
      hepmcvtx->add_particle_in(pp);
    }
  }

  /* Check that all particles (except beam particles) have a production vertex */
  for (unsigned p = 2; p < hepmc_particles.size(); p++)
  {
    if (!hepmc_particles.at(p)->production_vertex())
    {
      cout << "ReadEICFiles::process_event - Missing production vertex for one or more non-beam particles!" << endl;
      delete evt;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  /* Add HepMC vertices to event */
  for (unsigned v = 0; v < hepmc_vertices.size(); v++)
  {
    evt->add_vertex(hepmc_vertices.at(v));
  }

  /* set beam particles */
  evt->set_beam_particles(hepmc_beam1, hepmc_beam2);

  /* pass HepMC to PHNode*/
  PHHepMCGenEvent *success = PHHepMCGenHelper::insert_event(evt);
  if (Verbosity() > 1)
  {
    cout << __PRETTY_FUNCTION__ << " : " << __LINE__ << endl;
    evt->print();
    cout << endl << endl;
  }

  if (!success)
  {
    cout << "ReadEICFiles::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  /* Count up number of 'used' events from input file */
  entry++;

  /* Done */
  return 0;
}

int ReadEICFiles::CreateNodeTree(PHCompositeNode *topNode)
{
  PHHepMCGenHelper::create_node_tree(topNode);

  PHNodeIterator iter(topNode);
  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  EicEventHeader *evthead = findNode::getClass<EicEventHeader>(topNode, "EicEventHeader");
  if (!evthead)
  {
    evthead = new EicEventHeaderv1();
    PHIODataNode<PHObject> *HeaderNode = new PHIODataNode<PHObject>(evthead, "EicEventHeader", "PHObject");
    dstNode->addNode(HeaderNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
