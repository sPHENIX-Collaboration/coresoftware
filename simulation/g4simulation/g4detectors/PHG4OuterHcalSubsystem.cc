#include "PHG4OuterHcalSubsystem.h"

#include "PHG4HcalDefs.h"
#include "PHG4OuterHcalDetector.h"
#include "PHG4OuterHcalDisplayAction.h"
#include "PHG4OuterHcalSteppingAction.h"

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

#include <boost/foreach.hpp>

#include <cmath>     // for NAN
#include <iostream>  // for operator<<, basic_ostream
#include <set>       // for set
#include <sstream>

class PHG4Detector;

//_______________________________________________________________________
PHG4OuterHcalSubsystem::PHG4OuterHcalSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4OuterHcalSubsystem::~PHG4OuterHcalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4OuterHcalSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4OuterHcalDisplayAction(Name());

  // create detector
  m_Detector = new PHG4OuterHcalDetector(this, topNode, GetParams(), Name());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    std::set<std::string> nodes;
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }
    std::ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    nodes.insert(nodename.str());
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodename.str("");
      if (SuperDetector() != "NONE")
      {
        nodename << "G4HIT_ABSORBER_" << SuperDetector();
      }
      else
      {
        nodename << "G4HIT_ABSORBER_" << Name();
      }
      nodes.insert(nodename.str());
    }
    BOOST_FOREACH (std::string node, nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4OuterHcalSteppingAction(m_Detector, GetParams());
    m_SteppingAction->Init();
  }
  else
  {
    if (GetParams()->get_int_param("blackhole"))
    {
      m_SteppingAction = new PHG4OuterHcalSteppingAction(m_Detector, GetParams());
      m_SteppingAction->Init();
    }
  }

  return 0;
}

//_______________________________________________________________________
int PHG4OuterHcalSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4OuterHcalSubsystem::Print(const std::string &what) const
{
  std::cout << "Outer Hcal Parameters: " << std::endl;
  GetParams()->Print();
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4OuterHcalSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4OuterHcalSubsystem::SetLightCorrection(const double inner_radius, const double inner_corr, const double outer_radius, const double outer_corr)
{
  set_double_param("light_balance_inner_corr", inner_corr);
  set_double_param("light_balance_inner_radius", inner_radius);
  set_double_param("light_balance_outer_corr", outer_corr);
  set_double_param("light_balance_outer_radius", outer_radius);
  return;
}

void PHG4OuterHcalSubsystem::SetDefaultParameters()
{
  set_default_double_param("inner_radius", 183.3);
  set_default_double_param("light_balance_inner_corr", NAN);
  set_default_double_param("light_balance_inner_radius", NAN);
  set_default_double_param("light_balance_outer_corr", NAN);
  set_default_double_param("light_balance_outer_radius", NAN);
  // some math issue in the code does not subtract the magnet cutout correctly
  // (maybe some factor of 2 in a G4 volume creation)
  // The engineering drawing values are:
  //  set_default_double_param("magnet_cutout_radius", 195.31);
  //  set_default_double_param("magnet_cutout_scinti_radius", 195.96);
  // seting this to these values results in the correct edges
  // (verified by looking at the G4 hit coordinates of the inner edges)
  set_default_double_param("magnet_cutout_radius", 195.72);
  set_default_double_param("magnet_cutout_scinti_radius", 197.04);
  set_default_double_param("outer_radius", 264.71);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("scinti_eta_coverage", 1.1);
  set_default_double_param("scinti_gap", 0.85);
  set_default_double_param("scinti_gap_neighbor", 0.1);
  set_default_double_param("scinti_inner_radius", 183.89);
  // some math issue in the code subtracts 0.1mm+ so the scintillator
  // does not end at 263.27 as per drawing but at 263.26
  // adding 0.125mm compensates for this (so 263.2825 gives the desired 263.27
  set_default_double_param("scinti_outer_radius", 263.2825);
  set_default_double_param("scinti_tile_thickness", 0.7);
  set_default_double_param("size_z", 304.91 * 2);
  set_default_double_param("steplimits", NAN);
  set_default_double_param("tilt_angle", -11.23);  // engineering drawing
                                                   // corresponds very closely to 4 crossinge (-11.7826 deg)

  set_default_int_param("field_check", 0);
  set_default_int_param("light_scint_model", 1);
  set_default_int_param("magnet_cutout_first_scinti", 8);  // tile start at 0, drawing tile starts at 1

  // if ncross is set (and tilt_angle is NAN) tilt_angle is calculated
  // from number of crossings
  set_default_int_param("ncross", 0);
  set_default_int_param("n_towers", 64);
  set_default_int_param(PHG4HcalDefs::scipertwr, 5);
  set_default_int_param("n_scinti_tiles", 12);

  set_default_string_param("material", "Steel_1006");
  std::string defaultmapfilename;
  const char *Calibroot = getenv("CALIBRATIONROOT");
  if (Calibroot)
  {
    defaultmapfilename = Calibroot;
    defaultmapfilename += "/HCALOUT/tilemap/oHCALMaps092021.root";
  }
  set_default_string_param("MapFileName", defaultmapfilename);
  set_default_string_param("MapHistoName", "hCombinedMap");
}
