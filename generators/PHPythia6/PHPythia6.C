#include "PHPythia6.h"

#include <phhepmc/PHHepMCGenEvent.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <phool/PHIODataNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHObject.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHNodeReset.h>
#include <phool/PHTimeStamp.h>
#include <phool/PHRandomSeed.h>

#include <HepMC/PythiaWrapper.h>
#include <HepMC/IO_HEPEVT.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/IO_AsciiParticles.h>
#include <HepMC/GenEvent.h>

#include <gsl/gsl_randist.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <ctime>
#include <sys/time.h>
#include <algorithm>
#include <cctype>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;

PHPythia6::PHPythia6(const std::string &name):
  SubsysReco(name),
  _eventcount(0),
  _node_name("PHHepMCGenEvent"),
  _configFile("phpythia6.cfg"),
  _phhepmcevt(NULL),
  _save_ascii( false ),
  _filename_ascii("pythia_hepmc.dat"){

  //RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
}

PHPythia6::~PHPythia6() {
  //gsl_rng_free (RandomGenerator);
}

int PHPythia6::Init(PHCompositeNode *topNode) {

  /* Create node tree */
  CreateNodeTree(topNode);

  /* event numbering will start from 1 */
  _eventcount = 0;

  /* HepMC/example_MyPythia.cc:
   *
   * Pythia 6.1 uses HEPEVT with 4000 entries and 8-byte floating point
   *  numbers. We need to explicitly pass this information to the 
   *  HEPEVT_Wrapper.
   */
  HepMC::HEPEVT_Wrapper::set_max_number_entries(4000);
  HepMC::HEPEVT_Wrapper::set_sizeof_real(8);

  /* set pythia random number seed (mandatory!) */
  pydatr.mrpy[0]=55122 ;

  /* read pythia configuration and initialize */
  if (!_configFile.empty()) ReadConfig( _configFile );

  /* Initialization done */
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia6::End(PHCompositeNode *topNode) {

  //........................................TERMINATION
  // write out some information from Pythia to the screen
  call_pystat( 1 );

  return Fun4AllReturnCodes::EVENT_OK;
}

//__________________________________________________________
int PHPythia6::ReadConfig(const string cfg_file) {

  if ( cfg_file != "" ) _configFile = cfg_file;
  cout << "PHPythia6::read_config - Reading " << _configFile << endl;

  ifstream infile( _configFile );
  if (infile.fail ()) {
    cout << "PHPythia6::read_config - Failed to open file " << _configFile << endl;
    exit(2);
  }

  // initialize variables
  Int_t   _nevents(0);
  Float_t _roots(0);
  string  _proj;
  string  _targ;
  string  _frame;

  string FullLine;      // a complete line in the config file
  string label;         // the label

  int index = 999999;
  int ivalue = 999999;
  double value = 1e9;

  // get one line first
  getline(infile, FullLine);
  while ( !infile.eof() )
  {

    // skip lines that begin with #, or "//"
    if ( FullLine[0]=='#' || FullLine.substr(0, 2) == "//" )
    {
      getline(infile,FullLine);
      continue;
    }

    // make FullLine an istringstream
    istringstream line( FullLine.c_str() );

    // get label
    line >> label;

    // to lower case
    std::transform(label.begin(), label.end(), label.begin(), (int(*)(int)) std::tolower);

    // based on label, fill correct item
    if ( label == "nevents" )
    {
      line >> _nevents;
      cout << "nevents\t" << _nevents << endl;
    }
    else if ( label == "roots" )
    {
      line >> _roots;
      cout << "roots\t" << _roots << endl;
    }
    else if ( label == "proj" )
    {
      line >> _proj;
      cout << "proj\t" << _proj << endl;
    }
    else if ( label == "targ" )
    {
      line >> _targ;
      cout << "targ\t" << _targ << endl;
    }
    else if ( label == "frame" )
    {
      line >> _frame;
      cout << "frame\t" << _frame << endl;
    }
    else if ( label == "p1" || label == "p2")
    {
      int entry = 0;
      if ( label=="p1") entry = 1;
      else if ( label=="p2") entry = 2;

      char temp_line[10000];
      strcpy(temp_line,FullLine.c_str());
      char *sptr = strtok(temp_line," \t"); // skip first word
      sptr = strtok(NULL," \t");
      cout << label;
      int index = 1;
      while ( sptr != NULL )
        {
          double val = atof(sptr);
	  cout << "Entry: " << entry << endl;
	  index++;
          cout << "\t" << val;
          sptr = strtok(NULL," \t");
        }
      cout << endl;
    }
    else if ( label == "msel" )
      {
	line >> ivalue;
	pysubs.msel=ivalue;
	cout << "msel\t" << ivalue << endl;
      }
    else if ( label == "msub" )
    {
      line >> index >> ivalue;
      // careful with C/F77 differences: arrays in C start at 0, F77 at 1,
      // so we need to subtract 1 from the process #)
      pysubs.msub[index-1] = ivalue;
      cout << "msub\t" << index << " " << ivalue << endl;
    }
    else if ( label == "mstp" )
      {
	line >> index >> ivalue;
	pypars.mstp[index-1] = ivalue;
	cout << "mstp\t" << index << " " << ivalue << endl;
      }
    else if ( label == "mstj" )
      {
	line >> index >> ivalue;
	pydat1.mstj[index-1] = ivalue;
	cout << "mstj\t" << index << " " << ivalue << endl;
      }
    else if ( label == "mstu" )
      {
	line >> index >> ivalue;
	pydat1.mstu[index-1] = ivalue;
	cout << "mstu\t" << index << " " << ivalue << endl;
      }
    else if ( label == "ckin" )
      {
	line >> index >> value;
	pysubs.ckin[index-1] = ivalue;
	cout << "ckin\t" << index << " " << value << endl;
      }
    else if ( label == "parp" )
      {
	line >> index >> value;
	pypars.parp[index-1] = ivalue;
	cout << "parp\t" << index << " " << value << endl;
      }
    else if ( label == "parj" )
      {
	line >> index >> value;
	pydat1.parj[index-1] = ivalue;
	cout << "parj\t" << index << " " << value << endl;
      }
    else if ( label == "paru" )
      {
	line >> index >> value;
	pydat1.paru[index-1] = ivalue;
	cout << "paru\t" << index << " " << value << endl;
      }
    else if ( label == "parf" )
      {
	line >> index >> value;
	pydat2.parf[index-1] = ivalue;
	cout << "parf\t" << index << " " << value << endl;
      }
    else
      {
	// label was not understood
	cout << "************************************************************" << endl;
	cout << "PHPythia6::ReadConfig(), ERROR this option is not supported: " << FullLine << endl;
	cout << "************************************************************" << endl;
      }

    // get next line in file
    getline( infile, FullLine );
  }

  // Call pythia initialization
  call_pyinit( _frame.c_str(), _proj.c_str(), _targ.c_str(), _roots );

  infile.close();

  return _nevents;
}

//-* print pythia config info
void PHPythia6::print_config() const {
}

int PHPythia6::process_event(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "PHPythia6::process_event - event: " << _eventcount << endl;

  /* based on HepMC/example_MyPythia.cc
   *........................................HepMC INITIALIZATIONS
   *
   * Instantiate an IO strategy for reading from HEPEVT. */
  HepMC::IO_HEPEVT hepevtio;

  call_pyevnt();      // generate one event with Pythia
  // pythia pyhepc routine converts common PYJETS in common HEPEVT
  call_pyhepc( 1 );
  HepMC::GenEvent* evt = hepevtio.read_next_event();

  // define the units (Pythia uses GeV and mm)
  evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

  // add some information to the event
  evt->set_event_number(_eventcount);

  /* @TODO How to find out correct process ID from pythia? */
  //  evt->set_signal_process_id(20);

  // set number of multi parton interactions
  evt->set_mpi( pypars.msti[31-1] );

  // set cross section information
  evt->set_cross_section( HepMC::getPythiaCrossSection() );

  /* write the event out to the ascii files */
  if ( _save_ascii )
    {
      HepMC::IO_GenEvent ascii_io(_filename_ascii.c_str(),std::ios::out);
      ascii_io << evt;
    }

  /* pass HepMC to PHNode*/
  bool success = _phhepmcevt->addEvent(evt);
  if (!success) {
    cout << "PHPythia6::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /* print outs*/
  if (verbosity > 2) cout << "PHPythia6::process_event - FINISHED WHOLE EVENT" << endl;

  ++_eventcount;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia6::CreateNodeTree(PHCompositeNode *topNode) {

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

int PHPythia6::ResetEvent(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}
