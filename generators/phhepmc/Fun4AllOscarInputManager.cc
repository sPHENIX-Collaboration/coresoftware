#include "Fun4AllOscarInputManager.h"
#include "PHHepMCGenEventMap.h"

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/recoConsts.h>
#include <phool/getClass.h>

#include <ffaobjects/RunHeader.h>


#include <frog/FROG.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHObject.h>

#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>

#include <TString.h>
#include <TPRegexp.h>

#include <fstream>
#include <istream>
#include <iostream>
#include <sstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <cstdlib>
#include <memory>

using namespace std;

static boost::iostreams::filtering_streambuf<boost::iostreams::input> zinbuffer;
static const double toMM = 1.e-12;
typedef PHIODataNode<PHObject> PHObjectNode_t;

Fun4AllOscarInputManager::Fun4AllOscarInputManager(const string &name, const string &topnodename) :
  Fun4AllInputManager(name, ""),
  isopen(0),
  events_total(0),
  events_thisfile(0),
  topNodeName(topnodename),
  evt(NULL),
  skipEvents(0),
  skippedEvents(0),
  filestream(NULL),
  unzipstream(NULL),
  isCompressed(false)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  topNode = se->topNode(topNodeName.c_str());
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = se->getNode(InputNode(), topNodeName);

  PHHepMCGenEventMap *geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    geneventmap = new PHHepMCGenEventMap();
    PHIODataNode<PHObject> *newmapnode = new PHIODataNode<PHObject>(geneventmap, "PHHepMCGenEventMap", "PHObject");
    dstNode->addNode(newmapnode);
  }

  hepmc_helper.set_geneventmap(geneventmap);

}

Fun4AllOscarInputManager::~Fun4AllOscarInputManager()
{
  fileclose();
  delete filestream;
  delete unzipstream;
}

int
Fun4AllOscarInputManager::fileopen(const string &filenam)
{
  if (!mySyncManager)
    {
      cout << "Call fileopen only after you registered your Input Manager " << Name() << " with the Fun4AllServer" << endl;
      exit(1);
    }
  if (isopen)
    {
      cout << "Closing currently open file "
	   << filename
	   << " and opening " << filenam << endl;
      fileclose();
    }
  filename = filenam;
  FROG frog;
  string fname(frog.location(filename.c_str()));
  if (Verbosity() > 0)
    {
      cout << Name() << ": opening file " << fname << endl;
    }

  TString tstr(fname);
  TPRegexp bzip_ext(".bz2$");
  TPRegexp gzip_ext(".gz$");
  if (tstr.Contains(bzip_ext))
    {
      // use boost iosteam library to decompress bz2 on the fly
      filestream = new ifstream(fname.c_str(), std::ios::in | std::ios::binary);
      zinbuffer.push(boost::iostreams::bzip2_decompressor());
      zinbuffer.push(*filestream);
      unzipstream = new istream(&zinbuffer);
      isCompressed = true;
    }
  else if (tstr.Contains(gzip_ext))
    {
      // use boost iosream to decompress the gzip file on the fly
      filestream = new ifstream(fname.c_str(), std::ios::in | std::ios::binary);
      zinbuffer.push(boost::iostreams::gzip_decompressor());
      zinbuffer.push(*filestream);
      unzipstream = new istream(&zinbuffer);
      isCompressed = true;
    }
  else
    {
      theOscarFile.open(fname.c_str());
    }

  recoConsts *rc = recoConsts::instance();
  static bool run_number_forced = rc->FlagExist("RUNNUMBER");
  if ( run_number_forced )
    {
      mySyncManager->CurrentRun(rc->get_IntFlag("RUNNUMBER"));
    }
  else
    {
      mySyncManager->CurrentRun(-1);
    }
  events_thisfile = 0;
  isopen = 1;
  AddToFileOpened(fname); // add file to the list of files which were opened
  return 0;
}

int Fun4AllOscarInputManager::run(const int nevents)
{
 readagain:
  if (!isopen)
    {
      if (!filelist.size())
	{
	  if (Verbosity() > 0)
	    {
	      cout << Name() << ": No Input file open" << endl;
	    }
	  return 0;
	}
      else
	{
	  if (OpenNextFile())
	    {
	      cout << Name() << ": No Input file from filelist opened" << endl;
	      return 0;
	    }
	}
    }

  //Read oscar
  evt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
  int code = ConvertFromOscar();
  /*  
  if(skippedEvents < skipEvents || code == 3)
    {
      goto readagain;
    }
  */
  if (code == 1 || evt == NULL) 
    {
      if (Verbosity() > 1) cout << "Finished file!" << endl;
      fileclose();
      goto readagain;
    }

//  if(Verbosity() > 4) cout << "SIZE: " << phhepmcgenevt->size() << endl;
  //mySyncManager->CurrentEvent(evt->event_number());
  events_total++;
  events_thisfile++;

  // check if the local SubsysReco discards this event
  if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
    {
      //ResetEvent();
      goto readagain;
    }

  if(events_total < nevents) 
    {
      //ResetEvent();
      goto readagain;
    }

  return 0;
}

int
Fun4AllOscarInputManager::fileclose()
{
  if (!isopen)
    {
      cout << Name() << ": fileclose: No Input file open" << endl;
      return -1;
    }
  if(isCompressed)
    {
      filestream->close();
    }
  else
    {
      theOscarFile.close();
    }
  isopen = 0;
  // if we have a file list, move next entry to top of the list
  // or repeat the same entry again
  if (filelist.size() > 0)
    {
      if (repeat)
	{
	  filelist.push_back(*(filelist.begin()));
	  if (repeat > 0)
	    {
	      repeat--;
	    }
	}
      filelist.pop_front();
    }
  return 0;
}


void
Fun4AllOscarInputManager::Print(const string &what) const
{
  Fun4AllInputManager::Print(what);
  return ;
}

int
Fun4AllOscarInputManager::OpenNextFile()
{
  while (filelist.size() > 0)
    {
      list<string>::const_iterator iter = filelist.begin();
      if (Verbosity())
	{
	  cout << PHWHERE << " opening next file: " << *iter << endl;
	}
      if (fileopen((*iter).c_str()))
	{
	  cout << PHWHERE << " could not open file: " << *iter << endl;
	  filelist.pop_front();
	}
      else
	{
	  return 0;
	}

    }
  return -1;
}

int
Fun4AllOscarInputManager::ResetEvent()
{
  //delete evt;
  //evt = NULL;
  return 0;
}

int
Fun4AllOscarInputManager::PushBackEvents(const int i)
{
  int counter = 0;

  while(counter < i)
    {
      std::string theLine;
      while(getline(theOscarFile, theLine))
	{
	  if(theLine.find("#") == 0) continue;
	  vector<double> theInfo;
	  double number;
	  for(istringstream numbers_iss(theLine); numbers_iss >> number; )
	    {
	      theInfo.push_back(number);
	    }
	  
	  if(theInfo.size() == 2 && theInfo[0] == 0 && theInfo[1] == 0)
	    {
	      counter++;
	      skippedEvents++;
	      break;
	    }
	  else if (theInfo.size() == 2 && theInfo[0] == 0 && theInfo[1] > 0)
	    {
	      continue;
	    }
	}
      
    }

  if(theOscarFile.eof()) return -1;


  return 0;



}


int
Fun4AllOscarInputManager::ConvertFromOscar()
{

  delete evt;
  evt = NULL;
  
  if(theOscarFile.eof()) // if the file is exhausted bail out during this next read
    {
      cout << "Oscar EOF" << endl;
      return 1;
    }
  evt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);

  if(Verbosity() > 1) cout << "Reading Oscar Event " <<  events_total+skippedEvents+1 << endl;
  //Grab New Event From Oscar
  string theLine;
  vector< vector<double> > theEventVec;
  vector< HepMC::FourVector > theVtxVec;
  if(isCompressed)
    {
      // while(getline(unzipstream, theLine))
      // 	{
      // 	  if(theLine.find("#") == 0) continue;
      // 	  vector<double> theInfo; //format: N,pid,px,py,pz,E,mass,xvtx,yvtx,zvtx,?
      // 	  double number;
      // 	  for(istringstream numbers_iss(theLine); numbers_iss >> number; )
      // 	    {
      // 	      theInfo.push_back(number);
      // 	    }
	  
      // 	  if(theInfo.size() == 2 && theInfo[0] == 0 && theInfo[1] == 0)
      // 	    {
      // 	      break;
      // 	    }
      // 	  else if (theInfo.size() == 2 && theInfo[0] == 0 && theInfo[1] > 0)
      // 	    {
      // 	      continue;
      // 	    }
      // 	  else
      // 	    {
      // 	      theEventVec.push_back(theInfo);
      // 	      HepMC::FourVector vert(theInfo[8]*toMM, theInfo[9]*toMM, theInfo[10]*toMM, theInfo[11]);
      // 	      theVtxVec.push_back(vert);
      // 	    }
	  
      // 	}//while(getline)
    }
  else
    {
      while(getline(theOscarFile, theLine))
	{
	  if(theLine.find("#") == 0) continue;
	  vector<double> theInfo; //format: N,pid,px,py,pz,E,mass,xvtx,yvtx,zvtx,?
	  double number;
	  for(istringstream numbers_iss(theLine); numbers_iss >> number; )
	    {
	      theInfo.push_back(number);
	    }
	  
	  if(theInfo.size() == 2 && theInfo[0] == 0 && theInfo[1] == 0)
	    {
	      break;
	    }
	  else if (theInfo.size() == 2 && theInfo[0] == 0 && theInfo[1] > 0)
	    {
	      continue;
	    }
	  else
	    {
	      theEventVec.push_back(theInfo);
	    }
	  
	}//while(getline)

    }

  /*
  if(skippedEvents < skipEvents)
    {
      skippedEvents++;
      if (Verbosity() > 5) cout << "Skipping event " << skippedEvents << endl;
      return 2;
    }
  */

  //Set Event Number
  evt->set_event_number(events_total+1);

  //Loop Over One Event, Fill particles
  std::vector<HepMC::GenParticle*> hepevt_particles( theEventVec.size() );
  for(unsigned int i = 0; i < theEventVec.size(); i++)
    {
      //int N = (int)theEventVec[i][0];
      int pid = (int)theEventVec[i][1];
      double px = theEventVec[i][3];
      double py = theEventVec[i][4];
      double pz = theEventVec[i][5];
      double E = theEventVec[i][6];
      double m = theEventVec[i][7];
      int status = 1;//oscar only writes final state particles

      hepevt_particles[i] = new HepMC::GenParticle( HepMC::FourVector( px, py, pz, E ), pid, status );
      hepevt_particles[i]->setGeneratedMass(m);
      hepevt_particles[i]->suggest_barcode(i + 1);

    }

  for (unsigned int i = 0; i < theEventVec.size(); i++)
    {
      HepMC::GenParticle *p = hepevt_particles[i];
      HepMC::GenVertex* prod_vtx = p->production_vertex();
      if ( prod_vtx ) prod_vtx->add_particle_out( p );

      bool found = false;
      HepMC::FourVector prod_pos( theEventVec[i][8]*toMM, theEventVec[i][9]*toMM, theEventVec[i][10]*toMM, theEventVec[i][11] );
      if ( !prod_vtx )
        {
          //See if the vertex is already in the event 
          for(HepMC::GenEvent::vertex_iterator v = evt->vertices_begin(); v != evt->vertices_end(); ++v)
            {
	      HepMC::GenVertex *theV = *v;
              if(theV->position().x() != prod_pos.x()) continue;
              if(theV->position().y() != prod_pos.y()) continue;
              if(theV->position().z() != prod_pos.z()) continue;
              found = true;
              theV->add_particle_out(p);
            }
          if(!found)
            {
              //Didn't find vertex, add it
              prod_vtx = new HepMC::GenVertex(prod_pos);
              prod_vtx->add_particle_out( p );
              evt->add_vertex( prod_vtx );
            }
        }

      // If prod_vtx doesn't already have position specified, fill it.
      if ( !found && prod_vtx && prod_vtx->position() == HepMC::FourVector() ) prod_vtx->set_position( prod_pos );

    }

  
  evt->print();
  if(Verbosity() > 5) evt->print();
  if(Verbosity() > 3) cout << "Adding Event to phhepmcgenevt" << endl;

  PHHepMCGenEventMap::Iter ievt =
      hepmc_helper.get_geneventmap()->find(hepmc_helper.get_embedding_id());
  if (ievt != hepmc_helper.get_geneventmap()->end())
  {
    // override existing event
    ievt->second->addEvent(evt);
  }
  else
    hepmc_helper.insert_event(evt);
  return 0;

}

