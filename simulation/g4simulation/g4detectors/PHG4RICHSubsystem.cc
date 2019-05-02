// $$Id: PHG4RICHSubsystem.cc,v 1.1 2013/10/01 00:33:01 jinhuang Exp $$

/*!
 * \file PHG4RICHSubsystem.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.1 $$
 * \date $$Date: 2013/10/01 00:33:01 $$
 */

#include "PHG4RICHSubsystem.h"
#include "PHG4RICHDetector.h"
#include "PHG4RICHSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

using namespace ePHENIXRICH;

//_______________________________________________________________________
PHG4RICHSubsystem::PHG4RICHSubsystem( const char* name ):
PHG4Subsystem( name ),
detector_( NULL )
{
}

//_______________________________________________________________________
int PHG4RICHSubsystem::Init( PHCompositeNode* topNode )
{
  
  // create hit list
  PHG4HitContainer* rich_hits =  findNode::getClass<PHG4HitContainer>( topNode , "G4HIT_RICH" );
  if ( !rich_hits )
  {
    
    PHNodeIterator iter( topNode );
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST" ));
    dstNode->addNode( new PHIODataNode<PHObject>( rich_hits = new PHG4HitContainer("G4HIT_RICH"), "G4HIT_RICH","PHObject" ));
    
  }
  
  // create detector
  detector_ = new PHG4RICHDetector(topNode, geom);
  detector_->Verbosity(Verbosity());
  detector_->OverlapCheck(CheckOverlap());
  
  // create stepping action
  
  return 9;
  
}

//_______________________________________________________________________
int PHG4RICHSubsystem::process_event( PHCompositeNode* topNode )
{
  if ( PHG4RICHSteppingAction* p = dynamic_cast<PHG4RICHSteppingAction*>(detector_->GetSteppingAction()) )
    p->SetInterfacePointers( topNode );
  
  return 0;
  
}

//_______________________________________________________________________
PHG4Detector* PHG4RICHSubsystem::GetDetector( void ) const
{ return detector_; }

