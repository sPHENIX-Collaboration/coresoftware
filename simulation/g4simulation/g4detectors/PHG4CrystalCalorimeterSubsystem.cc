#include "PHG4CrystalCalorimeterSubsystem.h"
#include "PHG4CrystalCalorimeterDetector.h"
#include "PHG4CrystalCalorimeterDisplayAction.h"
#include "PHG4ProjCrystalCalorimeterDetector.h"
#include "PHG4CrystalCalorimeterSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;


//_______________________________________________________________________
PHG4CrystalCalorimeterSubsystem::PHG4CrystalCalorimeterSubsystem( const std::string &name, const int lyr ):
  PHG4Subsystem( name ),
  detector_( 0 ),
  m_SteppingAction( nullptr ),
  m_DisplayAction(nullptr),
  active(1),
  detector_type(name),
  mappingfile_(""),
  mappingfile_4x4_construct_(""),
  projective_(false)
{

}

//_______________________________________________________________________
PHG4CrystalCalorimeterSubsystem::~PHG4CrystalCalorimeterSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4CrystalCalorimeterSubsystem::Init( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create display settings before detector
  m_DisplayAction = new PHG4CrystalCalorimeterDisplayAction(Name());
  // create detector
  if ( projective_ )
    {
      cout << "PHG4CrystalCalorimeterSubsystem::InitRun - use PHG4ProjCrystalCalorimeterDetector" << endl;
      detector_ = new PHG4ProjCrystalCalorimeterDetector(topNode, Name());
      detector_->SetTowerMappingFile(mappingfile_);
      detector_->SetSupermoduleGeometry(mappingfile_4x4_construct_);
    }
  else
    {
      cout << "PHG4CrystalCalorimeterSubsystem::InitRun - use PHG4CrystalCalorimeterDetector" << endl;
      detector_ = new PHG4CrystalCalorimeterDetector(topNode, Name());
      detector_->SetTowerMappingFile(mappingfile_);
    }

  detector_->SetActive(active);
  detector_->SetAbsorberActive(active);
  detector_->OverlapCheck(CheckOverlap());

  if (active)
    {
      // create hit output node
      ostringstream nodename;
      nodename <<  "G4HIT_" << detector_type;

      PHG4HitContainer* crystal_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
      if (!crystal_hits)
        {
          crystal_hits = new PHG4HitContainer(nodename.str());
          PHIODataNode<PHObject> *hitNode = new PHIODataNode<PHObject>(crystal_hits, nodename.str().c_str(), "PHObject");
          dstNode->addNode(hitNode);
        }

      ostringstream absnodename;
      absnodename << "G4HIT_ABSORBER_" << detector_type;

      PHG4HitContainer* absorber_hits = findNode::getClass<PHG4HitContainer>(topNode, absnodename.str().c_str());
      if (!absorber_hits)
        {
          absorber_hits = new PHG4HitContainer(absnodename.str());
          PHIODataNode<PHObject> *abshitNode = new PHIODataNode<PHObject>(absorber_hits, absnodename.str().c_str(), "PHObject");
          dstNode->addNode(abshitNode);
        }

      // create stepping action
      m_SteppingAction = new PHG4CrystalCalorimeterSteppingAction(detector_);
    }
  return 0;
}


//_______________________________________________________________________
int
PHG4CrystalCalorimeterSubsystem::process_event( PHCompositeNode * topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
    {
      m_SteppingAction->SetInterfacePointers( topNode );
    }
  return 0;
}


//_______________________________________________________________________
PHG4Detector* PHG4CrystalCalorimeterSubsystem::GetDetector( void ) const
{
  return detector_;
}
