// $Id$

/*!
 * \file ${file_name}
 * \brief
 * \author Mickey Chiu <chiu@bnl.gov>
 * \version $Revision$
 * \date $Date$
 */

#include "PHG4PSTOFSubsystem.h"
#include "PHG4PSTOFDetector.h"
#include "PHG4ParametersContainer.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4PSTOFSteppingAction.h"

#include <g4main/PHG4HitContainer.h>

#include <pdbcalbase/PdbParameterMap.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <boost/foreach.hpp>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4PSTOFSubsystem::PHG4PSTOFSubsystem( const std::string &name, const int lyr ):
  PHG4DetectorGroupSubsystem( name, lyr ),
  detector_( NULL ),
  steppingAction_( NULL ),
  eventAction_(NULL)
{
  InitializeParameters();
}

//_______________________________________________________________________
int PHG4PSTOFSubsystem::InitRunSubsystem( PHCompositeNode* topNode )
{
  //cout << "In PHG4PSTOFSubsystem::InitRunSubsystem " << Name() << "\t" << GetParams()->get_int_param("active") << endl;
  cout << "In PHG4PSTOFSubsystem::InitRunSubsystem " << Name() << endl;
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4PSTOFDetector(topNode, GetParamsContainer(), Name());
  //detector_->SetActive(GetParamsContainer()->get_int_param("active"));
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());

  set<string> nodes;
  //if (GetParamsContainer()->get_int_param("active"))
  if ( 1 )  // need to figure out where the parameters are set
  {
    PHNodeIterator dstIter( dstNode );
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode",SuperDetector()));
    if (! DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }
    ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename <<  "G4HIT_" << SuperDetector();
      std::cout <<  "TOFXXXX G4HIT_" << Name() << endl;
    }
    else
    {
      nodename <<  "G4HIT_" << Name();
      //std::cout <<  "TOFYYYY G4HIT_" << Name() << endl;
    }
    nodes.insert(nodename.str());
    BOOST_FOREACH(string node, nodes)
    {
      PHG4HitContainer* g4_hits =  findNode::getClass<PHG4HitContainer>( topNode , node.c_str());
      if ( !g4_hits )
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode( new PHIODataNode<PHObject>( g4_hits, node.c_str(), "PHObject" ));
      }
      if (! eventAction_)
      {
        eventAction_ = new PHG4EventActionClearZeroEdep(topNode, node);
      }
      else
      {
        PHG4EventActionClearZeroEdep *evtact = dynamic_cast<PHG4EventActionClearZeroEdep *>(eventAction_);
        evtact->AddNode(node);
      }
    }
    // create stepping action
    steppingAction_ = new PHG4PSTOFSteppingAction(detector_, GetParamsContainer());
    steppingAction_->Init();
  }

  return 0;
}

//_______________________________________________________________________
int PHG4PSTOFSubsystem::process_event( PHCompositeNode * topNode )
{
  cout << "In PHG4PSTOFSubsystem::process_event()" << endl;

  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    cout << "Doing stepping action" << endl;

    steppingAction_->SetInterfacePointers( topNode );
  }
  return 0;
}

void PHG4PSTOFSubsystem::Print(const string &what) const
{
  //cout << "PSTOF Parameters: " << endl;
  PrintDefaultParams();
  PrintMacroParams();
  GetParamsContainer()->Print();
  if (detector_)
    {
      detector_->Print(what);
    }
  return;
}


//_______________________________________________________________________
PHG4Detector* PHG4PSTOFSubsystem::GetDetector( void ) const
{
    return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4PSTOFSubsystem::GetSteppingAction( void ) const
{
    return steppingAction_;
}

void PHG4PSTOFSubsystem::SetDefaultParameters()
{
  // whether to track through subsystem
  set_default_int_param(0,"active",1);

  // 
  set_default_int_param(0,"use_g4steps", 0);
 
  // geometry version number
  // we use negative numbers until the "official" version
  // when we build the detector
  set_default_int_param(0,"geometry_version",-1);

}

