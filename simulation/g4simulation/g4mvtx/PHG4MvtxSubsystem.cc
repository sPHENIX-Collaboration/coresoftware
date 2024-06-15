#include "PHG4MvtxSubsystem.h"

#include "PHG4MvtxDefs.h"
#include "PHG4MvtxDetector.h"
#include "PHG4MvtxDisplayAction.h"
#include "PHG4MvtxSteppingAction.h"

#include <mvtx/SegmentationAlpide.h>  // for Alpide constants

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4detectors/PHG4DetectorGroupSubsystem.h>  // for PHG4DetectorGrou...

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <cstdlib>   // for getenv
#include <iostream>  // for operator<<, basi...
#include <set>       // for _Rb_tree_const_i...
#include <sstream>
#include <utility>  // for pair

class PHG4Detector;

//_______________________________________________________________________
PHG4MvtxSubsystem::PHG4MvtxSubsystem(const std::string& name, const int _n_layers)
  : PHG4DetectorGroupSubsystem(name)
  , n_layers(_n_layers)
  , detector_type(name)
{
  for (unsigned int iLyr = 0; iLyr < n_layers; ++iLyr)
  {
    AddDetId(iLyr);
  }

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
  {
    std::cout << "PHG4MvtxSubsystem::Init started" << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector (detector adds its volumes to it)
  m_DisplayAction = new PHG4MvtxDisplayAction(Name());
  // create detector
  // These values are set from the calling macro using the setters defined in the .h file
  if (Verbosity())
  {
    std::cout << "    create Mvtx detector with " << n_layers << " layers." << std::endl;
  }
  m_Detector = new PHG4MvtxDetector(this, topNode, GetParamsContainer(), Name());
  m_Detector->Verbosity(Verbosity());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->Detector(detector_type);
  m_Detector->OverlapCheck(CheckOverlap());
  if (Verbosity())
  {
    std::cout << "    ------ created detector " << Name() << std::endl;
    GetParamsContainer()->Print();
  }
  // loop all layer to find atleast one active layer
  int active = 0;
  int supportactive = 0;
  int blackhole = 0;
  for (std::set<int>::const_iterator parContainerIter = GetDetIds().first; parContainerIter != GetDetIds().second; ++parContainerIter)
  {
    if (active || GetParamsContainer()->GetParameters(*parContainerIter)->get_int_param("active"))
    {
      active = 1;
    }
    if (supportactive || GetParamsContainer()->GetParameters(*parContainerIter)->get_int_param("supportactive"))
    {
      supportactive = 1;
    }
    if (blackhole || GetParamsContainer()->GetParameters(*parContainerIter)->get_int_param("blackhole"))
    {
      blackhole = 1;
    }
  }
  if (active)
  {
    std::set<std::string> nodes;
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* detNode = dstNode;
    if (SuperDetector() != "NONE" && !SuperDetector().empty())
    {
      detNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
      if (!detNode)
      {
        detNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(detNode);
      }
    }

    std::string detector_suffix = SuperDetector();
    if (detector_suffix == "NONE" || detector_suffix.empty())
    {
      detector_suffix = Name();
    }

    m_HitNodeName = "G4HIT_" + detector_suffix;
    nodes.insert(m_HitNodeName);
    m_SupportNodeName = "G4HIT_SUPPORT_" + detector_suffix;
    if (supportactive)
    {
      nodes.insert(m_SupportNodeName);
    }
    for (const auto& nodename : nodes)
    {
      PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename);
        if (Verbosity())
        {
          std::cout << PHWHERE << "creating hits node " << nodename << std::endl;
        }
        detNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename, "PHObject"));
      }
    }

    // create stepping action
    m_SteppingAction = new PHG4MvtxSteppingAction(m_Detector, GetParamsContainer());
    m_SteppingAction->SetHitNodeName("G4HIT", m_HitNodeName);
    m_SteppingAction->SetHitNodeName("G4HIT_SUPPORT", m_SupportNodeName);
    m_SteppingAction->Verbosity(Verbosity());
  }
  else
  {
    if (blackhole)
    {
      m_SteppingAction = new PHG4MvtxSteppingAction(m_Detector, GetParamsContainer());
    }
  }
  return 0;
}

//_______________________________________________________________________
int PHG4MvtxSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

//_______________________________________________________________________
PHG4Detector* PHG4MvtxSubsystem::GetDetector() const
{
  return m_Detector;
}

//_______________________________________________________________________
void PHG4MvtxSubsystem::SetDefaultParameters()
{
  for (std::set<int>::const_iterator lyr_it = GetDetIds().first; lyr_it != GetDetIds().second; ++lyr_it)
  {
    const int& ilyr = *lyr_it;
    const double rLr = PHG4MvtxDefs::mvtxdat[ilyr][PHG4MvtxDefs::kRmd];
    double turbo = radii2Turbo(PHG4MvtxDefs::mvtxdat[ilyr][PHG4MvtxDefs::kRmn], rLr, PHG4MvtxDefs::mvtxdat[ilyr][PHG4MvtxDefs::kRmx], SegmentationAlpide::SensorSizeRows * 10.);

    set_default_int_param(ilyr, "active", 1);  // non-automatic initialization in PHG4DetectorGroupSubsystem
    set_default_int_param(ilyr, "layer", ilyr);
    set_default_int_param(ilyr, "N_staves", PHG4MvtxDefs::mvtxdat[ilyr][PHG4MvtxDefs::kNStave]);

    set_default_double_param(ilyr, "layer_nominal_radius", rLr);
    set_default_double_param(ilyr, "phitilt", turbo);
    set_default_double_param(ilyr, "phi0", PHG4MvtxDefs::mvtxdat[ilyr][PHG4MvtxDefs::kPhi0]);
    set_default_string_param(ilyr, "material", "G4_AIR");  // default - almost nothing
  }

  set_default_string_param(PHG4MvtxDefs::GLOBAL, "stave_geometry_file", "ITS.gdml");  // default - almost nothing
  char* calibrationsroot = getenv("CALIBRATIONROOT");
  std::string end_wheels_sideS = "ITS_ibEndWheelSideA.gdml";
  std::string end_wheels_sideN = "ITS_ibEndWheelSideC.gdml";
  if (calibrationsroot != nullptr)
  {
    end_wheels_sideS = std::string(calibrationsroot) + std::string("/Tracking/geometry/") + end_wheels_sideS;
    end_wheels_sideN = std::string(calibrationsroot) + std::string("/Tracking/geometry/") + end_wheels_sideN;
  }
}
