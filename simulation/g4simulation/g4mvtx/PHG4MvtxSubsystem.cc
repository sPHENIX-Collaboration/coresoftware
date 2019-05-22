#include "PHG4MvtxSubsystem.h"

#include "PHG4MvtxDefs.h"
#include "PHG4MvtxDetector.h"
#include "PHG4MvtxDisplayAction.h"
#include "PHG4MvtxSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <Geant4/G4GDMLParser.hh>

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4MvtxSubsystem::PHG4MvtxSubsystem(const std::string& name, const int _n_layers)
  : PHG4DetectorGroupSubsystem(name)
  , m_Detector(nullptr)
  , steppingAction_(nullptr)
  , m_DisplayAction(nullptr)
  , n_layers(_n_layers)
  , detector_type(name)
{
  for (unsigned int iLyr = 0; iLyr < n_layers; ++iLyr)
    AddDetId(iLyr);

  InitializeParameters();

  // put the layers into name so we get unique names
  // for multiple layers
  Name(name);
  SuperDetector(name);
}

//_______________________________________________________________________
PHG4MvtxSubsystem::~PHG4MvtxSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4MvtxSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
    cout << "PHG4MvtxSubsystem::Init started" << endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector (detector adds its volumes to it)
  m_DisplayAction = new PHG4MvtxDisplayAction(Name());
  // create detector
  // These values are set from the calling macro using the setters defined in the .h file
  if (Verbosity())
  {
    cout << "    create Mvtx detector with " << n_layers << " layers." << endl;
  }
  m_Detector = new PHG4MvtxDetector(this, topNode, GetParamsContainer(), Name());
  m_Detector->Verbosity(Verbosity());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->Detector(detector_type);
  m_Detector->OverlapCheck(CheckOverlap());
  if (Verbosity())
  {
    cout << "    ------ created detector " << Name() << endl;
    GetParamsContainer()->Print();
  }
  //loop all layer to find atleast one active layer
  int active = 0;
  // for now not absorber are implemnented yet
  int absorberactive = 0;
  int blackhole = 0;
  for (set<int>::const_iterator parContainerIter = GetDetIds().first; parContainerIter != GetDetIds().second; ++parContainerIter)
  {
    if (active || GetParamsContainer()->GetParameters(*parContainerIter)->get_int_param("active"))
    {
      active = 1;
    }
    if (absorberactive || GetParamsContainer()->GetParameters(*parContainerIter)->get_int_param("absorberactive"))
    {
      absorberactive = 1;
    }
    if (blackhole || GetParamsContainer()->GetParameters(*parContainerIter)->get_int_param("blackhole"))
    {
      blackhole = 1;
    }
  }
  if (active)
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* detNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!detNode)
    {
      detNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(detNode);
    }
    ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << detector_type;
    }
    // create hit list
    PHG4HitContainer* block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!block_hits)
    {
      detNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }
    if (Verbosity())
      cout << PHWHERE << "creating hits node " << nodename.str() << endl;

    if (absorberactive)
    {
      nodename.str("");
      if (SuperDetector() != "NONE")
      {
        nodename << "G4HIT_ABSORBER_" << SuperDetector();
      }
      else
      {
        nodename << "G4HIT_ABSORBER_" << detector_type;
      }
      block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
      if (!block_hits)
      {
        detNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
      }
      if (Verbosity())
        cout << PHWHERE << "creating hits node " << nodename.str() << endl;
    }

    // create stepping action
    steppingAction_ = new PHG4MvtxSteppingAction(m_Detector);
  }
  else
  {
    if (blackhole)
    {
      steppingAction_ = new PHG4MvtxSteppingAction(m_Detector);
    }
  }
  return 0;
}

//_______________________________________________________________________
int PHG4MvtxSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers(topNode);
  }
  return 0;
}

//_______________________________________________________________________
PHG4Detector* PHG4MvtxSubsystem::GetDetector(void) const
{
  return m_Detector;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4MvtxSubsystem::GetSteppingAction(void) const
{
  return steppingAction_;
}

void PHG4MvtxSubsystem::SetDefaultParameters()
{
  //TODO: Move to defMvtx at some point
  const int kNLr = 3;
  enum
  {
    kRmn,
    kRmd,
    kRmx,
    kNModPerStave,
    kPhi0,
    kNStave,
    kNPar
  };
  // Radii are from last TDR (ALICE-TDR-017.pdf Tab. 1.1, rMid is mean value)
  const double mvtxdat[kNLr][kNPar] = {
      {22.4, 23.4, 26.7, 9., 0., 12.},  // for each layer: rMin,rMid,rMax,NChip/Stave, phi0, nStaves
      {30.1, 31.5, 34.6, 9., 0., 16.},
      {37.8, 39.3, 42.1, 9., 0., 20.}};

  for (set<int>::const_iterator lyr_it = GetDetIds().first; lyr_it != GetDetIds().second; ++lyr_it)
  {
    const int& ilyr = *lyr_it;
    set_default_int_param(ilyr, "active", 1);  //non-automatic initialization in PHG4DetectorGroupSubsystem
    set_default_int_param(ilyr, "layer", ilyr);
    set_default_int_param(ilyr, "N_staves", mvtxdat[ilyr][kNStave]);

    set_default_double_param(ilyr, "layer_nominal_radius", mvtxdat[ilyr][kRmd]);
    set_default_double_param(ilyr, "phitilt", NAN);

    set_default_string_param(ilyr, "material", "G4_AIR");  // default - almost nothing
  }
  set_default_string_param(PHG4MvtxDefs::GLOBAL, "stave_geometry_file", "ITS.gdml");  // default - almost nothing
  set_default_double_param(PHG4MvtxDefs::ALPIDE_SEGMENTATION, "pixel_x", NAN);
  set_default_double_param(PHG4MvtxDefs::ALPIDE_SEGMENTATION, "pixel_z", NAN);
  set_default_double_param(PHG4MvtxDefs::ALPIDE_SEGMENTATION, "pixel_thickness", NAN);
}
