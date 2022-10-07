#include "PHG4IHCalSubsystem.h"

#include "PHG4IHCalDetector.h"
#include "PHG4IHCalDisplayAction.h"
#include "PHG4IHCalSteppingAction.h"

#include <g4detectors/PHG4DetectorSubsystem.h>  // for PHG4DetectorSubsystem
#include <g4detectors/PHG4HcalDefs.h>

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>


#include <cmath>     // for NAN
#include <iostream>  // for operator<<, basic_ostream
#include <set>

class PHG4Detector;

//_______________________________________________________________________
PHG4IHCalSubsystem::PHG4IHCalSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4IHCalSubsystem::~PHG4IHCalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4IHCalSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4IHCalDisplayAction(Name());

  // create detector
  m_Detector = new PHG4IHCalDetector(this, topNode, GetParams(), Name());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  std::set<std::string> nodes;
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dstNode;
   
    if (SuperDetector() != "NONE" && !SuperDetector().empty())
    {
      PHNodeIterator iter_dst(dstNode);
      DetNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));

      if (!DetNode)
      {
        DetNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(DetNode);
      }
    }
    std::string detector_suffix = SuperDetector();
    if (detector_suffix == "NONE" || detector_suffix.empty())
    {
      detector_suffix = Name();
    }

    m_HitNodeName = "G4HIT_" + detector_suffix;
    nodes.insert(m_HitNodeName);
    m_AbsorberNodeName = "G4HIT_ABSORBER_" + detector_suffix;
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodes.insert(m_AbsorberNodeName);
    }
    for (const auto &nodename : nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename, "PHObject"));
      }
    }

    // create stepping action
    m_SteppingAction = new PHG4IHCalSteppingAction(m_Detector, GetParams());
    m_SteppingAction->Init();
    m_SteppingAction->SetHitNodeName("G4HIT", m_HitNodeName);
    m_SteppingAction->SetHitNodeName("G4HIT_ABSORBER", m_AbsorberNodeName);
  }
  else
  {
    // if this is a black hole it does not have to be active
    if (GetParams()->get_int_param("blackhole"))
    {
      m_SteppingAction = new PHG4IHCalSteppingAction(m_Detector, GetParams());
      m_SteppingAction->Init();
    }
  }
  return 0;
}

//_______________________________________________________________________
int PHG4IHCalSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4IHCalSubsystem::Print(const std::string &what) const
{
  std::cout << Name() << " Parameters: " << std::endl;
  GetParams()->Print();
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
  if (m_SteppingAction)
  {
    m_SteppingAction->Print(what);
  }

  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4IHCalSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4IHCalSubsystem::SetDefaultParameters()
{
  set_default_double_param(PHG4HcalDefs::innerrad, 115);
  set_default_double_param("light_balance_inner_corr", NAN);
  set_default_double_param("light_balance_inner_radius", NAN);
  set_default_double_param("light_balance_outer_corr", NAN);
  set_default_double_param("light_balance_outer_radius", NAN);
  set_default_double_param(PHG4HcalDefs::outerrad, 274.010 / 2 + 3);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 180.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("size_z", 435.000 + 10 );
  set_default_double_param("Birk_const", 0.07943);

  set_default_int_param("light_scint_model", 1);
  set_default_int_param(PHG4HcalDefs::n_towers, 64);
  set_default_int_param(PHG4HcalDefs::scipertwr, 4);
  set_default_int_param(PHG4HcalDefs::n_scinti_tiles, 12);

  set_default_string_param("GDMPath", "DefaultParameters-InvadPath");
  const char* Calibroot = getenv("CALIBRATIONROOT");
  std::string defaultmapfilename;
  if (Calibroot)
  {
    defaultmapfilename = Calibroot;
    defaultmapfilename += "/HCALIN/tilemap/ihcalgdmlmap09212022.root";
  }
  set_default_string_param("MapFileName", defaultmapfilename);
  set_default_string_param("MapHistoName", "ihcalcombinedgdmlnormtbyt");
}

void PHG4IHCalSubsystem::SetLightCorrection(const double inner_radius, const double inner_corr, const double outer_radius, const double outer_corr)
{
  set_double_param("light_balance_inner_corr", inner_corr);
  set_double_param("light_balance_inner_radius", inner_radius);
  set_double_param("light_balance_outer_corr", outer_corr);
  set_double_param("light_balance_outer_radius", outer_radius);
  return;
}
