#include "PHPythia6.h"
#include "PHPy6GenTrigger.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

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

#define pytune pytune_
extern "C" int  pytune(int *itune);

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;

PHPythia6::PHPythia6(const std::string &name):
  SubsysReco(name),
  _eventcount(0),
  _geneventcount(0),
  _configFile("phpythia6.cfg"),
  _save_ascii( false ),
  _filename_ascii("pythia_hepmc.dat"),
  _registeredTriggers(),
  _triggersOR(true),
  _triggersAND(false){

  hepmc_helper.set_embedding_id(1); // default embedding ID to 1
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
  int fSeed = PHRandomSeed(); 
  fSeed = abs(fSeed)%900000000;
	    	  	  
  if ( (fSeed>=0) && (fSeed<=900000000) ) {
    pydatr.mrpy[0] = fSeed;                   // set seed
    pydatr.mrpy[1] = 0;                       // use new seed
  } else {
    cout << PHWHERE << " ERROR: seed " << fSeed << " is not valid" << endl;
    exit(2); 
  }

  /* read pythia configuration and initialize */
  if (!_configFile.empty()) ReadConfig( _configFile );

  /* Initialization done */
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia6::End(PHCompositeNode *topNode) {

  //........................................TERMINATION
  // write out some information from Pythia to the screen
  call_pystat( 1 );

  //match pythia printout
  cout << " |                                                                "
       << "                                                 | " << endl; 
  cout << "                         PHPythia6::End - " << _eventcount
       << " events passed trigger" << endl;
  cout << "                         Fraction passed: " << _eventcount
       << "/" << _geneventcount
       << " = " << _eventcount/float(_geneventcount) << endl;
  cout << " *-------  End PYTHIA Trigger Statistics  ------------------------"
       << "-------------------------------------------------* " << endl;

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
      if (label == "p1") //Momentum of Projectile Beam (e- for e-p)
	{
	  line >> index >> value;
	  //Index Options(3MOM): 1 = x-momentum; 2 = y-momentum; 3 = z-momentum
	  pyjets.p[index-1][0] = value;
	  cout << "p1\t" << index << " " << value << endl;
	}
      if (label == "p2") //Momentum of Target Beam (p for e-p)
	{
	  line >> index >> value;
	  //Index Options(3MOM): 1 = x-momentum; 2 = y-momentum; 3 = z-momentum
	  pyjets.p[index-1][1] = value;
	  cout << "p2\t" << index << " " << value << endl;
	}
      /*
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
	cout << endl;*/
    }
    else if ( label == "msel" )
      {
	line >> value;
	pysubs.msel= (int) value;
	cout << "msel\t" << value << endl;
	IntegerTest(value);
      }
    else if ( label == "msub" )
    {
      line >> index >> value;
      // careful with C/F77 differences: arrays in C start at 0, F77 at 1,
      // so we need to subtract 1 from the process #)
      pysubs.msub[index-1] = (int) value;
      cout << "msub\t" << index << " " << value << endl;
      IntegerTest(value);
    }
    else if ( label == "mstp" )
      {
	line >> index >> value;
	pypars.mstp[index-1] = (int) value;
	cout << "mstp\t" << index << " " << value << endl;
	IntegerTest(value);
      }
    else if ( label == "mstj" )
      {
	line >> index >> value;
	pydat1.mstj[index-1] = (int) value;
	cout << "mstj\t" << index << " " << value << endl;
	IntegerTest(value);
      }
    else if ( label == "mstu" )
      {
	line >> index >> value;
	pydat1.mstu[index-1] = (int) value;
	cout << "mstu\t" << index << " " << value << endl;
	IntegerTest(value);
      }
    else if ( label == "ckin" )
      {
	line >> index >> value;
	pysubs.ckin[index-1] =value;
	cout << "ckin\t" << index << " " << value << endl;
      }
    else if ( label == "parp" )
      {
	line >> index >> value;
	pypars.parp[index-1] = value;
	cout << "parp\t" << index << " " << value << endl;
      }
    else if ( label == "parj" )
      {
	line >> index >> value;
	pydat1.parj[index-1] = value;
	cout << "parj\t" << index << " " << value << endl;
      }
    else if ( label == "paru" )
      {
	line >> index >> value;
	pydat1.paru[index-1] = value;
	cout << "paru\t" << index << " " << value << endl;
      }
    else if ( label == "parf" )
      {
	line >> index >> value;
	pydat2.parf[index-1] = value;
	cout << "parf\t" << index << " " << value << endl;
      }
    else if ( label == "mdme" )
    {
      int idc = 0;          // decay channel
      line >> idc >> index >> value;
      
      // if (ivalue==1/0) turn on/off decay channel idc
      pydat3.mdme[index-1][idc-1] = value; 
      cout << "mdme\t" << idc << " " << index << " " << value << endl;
    }
    else if ( label == "pmas" )
    {
      int idc = 0;          
      line >> idc >> index >> value;
      
      pydat2.pmas[index-1][idc-1] = value; 
      cout << "pmas\t" << idc << " " << index << " " << value << endl;
    }
    else if ( label == "pytune" )
    {
      int ivalue; 
      line >> ivalue;
      pytune(&ivalue);
      cout << "pytune\t" << ivalue << endl;
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

  //call_pylist(12); 

  infile.close();

  return _nevents;
}

//-* print pythia config info
void PHPythia6::print_config() const {
}

int PHPythia6::process_event(PHCompositeNode *topNode) {

  if (Verbosity() > 1) cout << "PHPythia6::process_event - event: " << _eventcount << endl;

  bool passedTrigger = false;
  std::vector<bool> theTriggerResults;
  int genCounter = 0;

  /* based on HepMC/example_MyPythia.cc
   *........................................HepMC INITIALIZATIONS
   *
   * Instantiate an IO strategy for reading from HEPEVT. */
  HepMC::IO_HEPEVT hepevtio;
  HepMC::GenEvent* evt; 

  while (!passedTrigger) {
    ++genCounter;

    call_pyevnt();      // generate one event with Pythia
    _geneventcount++; 
    // pythia pyhepc routine converts common PYJETS in common HEPEVT
    call_pyhepc( 1 );
    evt = hepevtio.read_next_event();

    // define the units (Pythia uses GeV and mm)
    evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

    // add some information to the event
    evt->set_event_number(_eventcount);

    /* process ID from pythia */
    evt->set_signal_process_id(pypars.msti[1-1]);

    // set number of multi parton interactions
    evt->set_mpi( pypars.msti[31-1] );

    // set cross section information
    evt->set_cross_section( HepMC::getPythiaCrossSection() );

    // Set the PDF information
    HepMC::PdfInfo pdfinfo;
    pdfinfo.set_x1(pypars.pari[33-1]); 
    pdfinfo.set_x2(pypars.pari[34-1]); 
    pdfinfo.set_scalePDF(pypars.pari[22-1]);
    pdfinfo.set_id1(pypars.msti[15-1]); 
    pdfinfo.set_id2(pypars.msti[16-1]);  
    evt->set_pdf_info(pdfinfo); 

    // test trigger logic
    
    bool andScoreKeeper = true;
    if (Verbosity() > 2) {
      cout << "PHPythia6::process_event - triggersize: " << _registeredTriggers.size() << endl;
    }

    for (unsigned int tr = 0; tr < _registeredTriggers.size(); tr++) { 
      bool trigResult = _registeredTriggers[tr]->Apply(evt);

      if (Verbosity() > 2) {
	cout << "PHPythia6::process_event trigger: "
	     << _registeredTriggers[tr]->GetName() << "  " << trigResult << endl;
      }

      if (_triggersOR && trigResult) {
	passedTrigger = true;
	break;
      } else if (_triggersAND) {
	andScoreKeeper &= trigResult;
      }
      
      if (Verbosity() > 2 && !passedTrigger) {
	cout << "PHPythia8::process_event - failed trigger: "
	     << _registeredTriggers[tr]->GetName() <<  endl;
      }
    }

    if ((andScoreKeeper && _triggersAND) || (_registeredTriggers.size() == 0)) {
      passedTrigger = true;
      genCounter = 0;
    }
    
    // delete failed events
    if(!passedTrigger) delete evt; 

  }

  /* write the event out to the ascii files */
  if ( _save_ascii )
  {
    HepMC::IO_GenEvent ascii_io(_filename_ascii.c_str(),std::ios::out);
    ascii_io << evt;
  }

  /* pass HepMC to PHNode*/

  PHHepMCGenEvent * success = hepmc_helper . insert_event(evt);
  if (!success) {
    cout << "PHPythia6::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  /* print outs*/
  if (Verbosity() > 2) cout << "PHPythia6::process_event - FINISHED WHOLE EVENT" << endl;

  ++_eventcount;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia6::CreateNodeTree(PHCompositeNode *topNode) {

  hepmc_helper.create_node_tree(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHPythia6::IntegerTest(double number ) {

  if (fmod(number, 1.0) != 0) {
    cout << "Warning: Value " << number << " is not an integer." << endl;
    cout << "This parameter requires an integer value." << endl;
    cout << "Value of parameter truncated to " << (int) number  << endl;

    //...End simulation if a double value is input for an integer parameter
    //    throw Fun4AllReturnCodes::ABORTRUN;
  }
  return;
}

int PHPythia6::ResetEvent(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHPythia6::register_trigger(PHPy6GenTrigger *theTrigger) {
  if(Verbosity() > 1) cout << "PHPythia6::registerTrigger - trigger " << theTrigger->GetName() << " registered" << endl;
  _registeredTriggers.push_back(theTrigger);
}
