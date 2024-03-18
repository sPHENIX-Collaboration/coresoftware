#include "XingShiftCal.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/recoConsts.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/msg_profile.h>

#include <TH1.h>
#include <TCanvas.h>

#include <phool/PHCompositeNode.h>


XingShiftCal::XingShiftCal(const std::string &name):
 SubsysReco(name)
{

  nevt = 0;
  for (int i = 0; i < NTRIG; i++){
    for (int j = 0; j < NBUNCHES; j++){
      scalercounts[i][j] = 0;
    }
  }
  std::cout << "XingShiftCal::XingShiftCal(const std::string &name) Calling ctor" << std::endl;
}


XingShiftCal::~XingShiftCal()
{
  std::cout << "XingShiftCal::~XingShiftCal() Calling dtor" << std::endl;
}


int XingShiftCal::Init(PHCompositeNode *topNode)
{
  std::cout << "XingShiftCal::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


int XingShiftCal::InitRun(PHCompositeNode *topNode)
{

  //std::cout << "XingShiftCal::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


int XingShiftCal::process_event(PHCompositeNode *topNode)
{


  if (done == 1){
    return Fun4AllReturnCodes::EVENT_OK;
  }
  //std::cout << "XingShiftCal::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  Event *evt = findNode::getClass<Event>(topNode,"PRDF");
  
  /*
  if (evt->getEvtType() == BEGRUNEVENT){
    pBlueSpin = evt->getPacket(packet_BLUESPIN);
    pYellSpin = evt->getPacket(packet_YELLSPIN);
    for (int i = 0; i < NBUNCHES; i++){
      blueSpinPattern[i] = pBlueSpin->iValue(i);
      yellSpinPattern[i] = pYellSpin->iValue(i);
    }
    delete pBlueSpin;
    delete pYellSpin;
    //get_cdev_info(p);
    
  } else if 
  */
  if (evt->getEvtType() == DATAEVENT){
    p = evt->getPacket(packet_GL1);
    int bunchnr = p->lValue(0,"BunchNumber");
    for (int i = 0; i < NTRIG; i++){
      //2nd arg of lValue: 0 is raw trigger count, 1 is live trigger count, 2 is scaled trigger count
      // this is the gl1 scaler for now until gl1p scaler is implemented in the GL1 packets
      long gl1pscaler = p->lValue(i,2);
      scalercounts[i][bunchnr] = gl1pscaler;      
    }
    delete p;
  }

  
  if (nevt > threshold){
    Calibrate();
    if (success){
      done = 1;
    } else {
      threshold += threshold;
    }
  }

  if (nevt > evtcap){
    done = 1;
  }

  nevt++;
  
  return Fun4AllReturnCodes::EVENT_OK;
}


int XingShiftCal::ResetEvent(PHCompositeNode *topNode)
{
  //std::cout << "XingShiftCal::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


int XingShiftCal::EndRun(const int runnumber)
{
  std::cout << "XingShiftCal::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


int XingShiftCal::End(PHCompositeNode *topNode)
{
  std::cout << "XingShiftCal::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  if (done){
    //if (!success){
      Calibrate(1);
      //}
    //CommitConstantsToFile();
  } 
  else {
    success = 0;
    std::cout << "Not enough statistics. Did not calibrate." << std::endl;
  }
  
  /*
  if (committopdbcal && success){
    CommitConstantsToDatabase();
    CommitPatternToSpinDB();
  }
  */

  return Fun4AllReturnCodes::EVENT_OK;
}


int XingShiftCal::Calibrate(const int final)
{  
  CalculateCrossingShift(xingshift,scalercounts,success);
  if (!success){
    std::cout << "CROSSING CALIBRATION FAILED." << std::endl;
    if (final){std::cout << "DONE CALIBRATING." << std::endl;}
    return 0;
  }
  
  if (final)
    {
      //ostringstream comment;
      std::cout << "CROSSING CALIBRATION SUCCESS. XINGSHIFT: " << xingshift << std::endl;
      std::cout << "DONE CALIBRATING." << std::endl;
      //AddComment(comment.str());
    }

  
  
  return 0;
}

int XingShiftCal::CalculateCrossingShift(int& xingshift, long long counts[NTRIG][NBUNCHES], bool& success){
  
  success = false;
  int shift_array[NTRIG] = {0};

  int trig_inactive_array[NTRIG] = {0};
  
  int _temp;
  for(int itrig = 0; itrig < NTRIG; itrig++){
    
    long long _counts = 0;
    for(int ii = 0; ii< NBUNCHES; ii++){
      _counts += counts[itrig][ii];
    }
    
    if (_counts < 10000){trig_inactive_array[itrig] = 1;}

    long long abort_sum_prev = _counts;

    _temp = 0;
    for (int ishift = 0; ishift < NBUNCHES; ishift++){
      long long abort_sum = 0;
      for (int iabortbunch = NBUNCHES-9; iabortbunch < NBUNCHES; iabortbunch++){ 
	abort_sum += counts[itrig][(iabortbunch+ishift)%NBUNCHES];
      }
      if (abort_sum < abort_sum_prev){
	abort_sum_prev = abort_sum;
	_temp = ishift;
      }
    }
  
    shift_array[itrig] = _temp;

    if ( itrig > 0) { //if not matching for all trigger selections used, fails
      if ( shift_array[itrig]!= shift_array[itrig-1] && !trig_inactive_array[itrig] && !trig_inactive_array[itrig-1]) { 
	xingshift = 0;
	success = false;
	return 0;
      } else if (!trig_inactive_array[itrig]){
	xingshift = shift_array[itrig];
      }
    }
  }

  success = true;
  return 0;
}



//____________________________________________________________________________..
int XingShiftCal::Reset(PHCompositeNode *topNode)
{
  //std::cout << "XingShiftCal::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void XingShiftCal::Print(const std::string &what) const
{
  std::cout << "XingShiftCal::Print(const std::string &what) const Printing info for " << what << std::endl;
}
