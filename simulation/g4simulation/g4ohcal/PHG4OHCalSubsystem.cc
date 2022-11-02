#include "PHG4OHCalSubsystem.h"

#include "PHG4OHCalDetector.h"
#include "PHG4OHCalDisplayAction.h"
#include "PHG4OHCalSteppingAction.h"

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
#include <set>       // for set
#include <cstdlib> // for getenv
class PHG4Detector;

//_______________________________________________________________________
PHG4OHCalSubsystem::PHG4OHCalSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4OHCalSubsystem::~PHG4OHCalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4OHCalSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4OHCalDisplayAction(Name());

  if (get_string_param("IronFieldMapPath") == "DefaultParameters-InvadPath" )
  {
    std::cout <<__PRETTY_FUNCTION__<<": invalid string parameter IronFieldMapPath, where we expect a 3D field map"<<std::endl;
    exit (1);
  }


  // create detector
  m_Detector = new PHG4OHCalDetector(this, topNode, GetParams(), Name());
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
    m_SteppingAction = new PHG4OHCalSteppingAction(m_Detector, GetParams());
    m_SteppingAction->Init();
    m_SteppingAction->SetHitNodeName("G4HIT", m_HitNodeName);
    m_SteppingAction->SetHitNodeName("G4HIT_ABSORBER", m_AbsorberNodeName);
  }
  else
  {
    if (GetParams()->get_int_param("blackhole"))
    {
      m_SteppingAction = new PHG4OHCalSteppingAction(m_Detector, GetParams());
      m_SteppingAction->Init();
    }
  }

  return 0;
}

//_______________________________________________________________________
int PHG4OHCalSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4OHCalSubsystem::Print(const std::string &what) const
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
PHG4Detector *PHG4OHCalSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4OHCalSubsystem::SetLightCorrection(const double inner_radius, const double inner_corr, const double outer_radius, const double outer_corr)
{
  set_double_param("light_balance_inner_corr", inner_corr);
  set_double_param("light_balance_inner_radius", inner_radius);
  set_double_param("light_balance_outer_corr", outer_corr);
  set_double_param("light_balance_outer_radius", outer_radius);
  return;
}

void PHG4OHCalSubsystem::SetDefaultParameters()
{
  set_default_double_param("inner_radius", 182.423  - 5);
  set_default_double_param("light_balance_inner_corr", NAN);
  set_default_double_param("light_balance_inner_radius", NAN);
  set_default_double_param("light_balance_outer_corr", NAN);
  set_default_double_param("light_balance_outer_radius", NAN);
  set_default_double_param("outer_radius", 269.317  + 5 );
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 180.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("size_z", 639.240 + 10);
  set_default_double_param("Birk_const", 0.07943);
  set_default_int_param("field_check", 0);
  set_default_int_param("light_scint_model", 1);

  set_default_int_param("n_towers", 64);
  set_default_int_param(PHG4HcalDefs::scipertwr, 5);
  set_default_int_param("n_scinti_tiles", 12);

  set_default_string_param("GDMPath", "DefaultParameters-InvadPath");
  std::string defaultmapfilename;
  const char* Calibroot = getenv("CALIBRATIONROOT");
  if (Calibroot)
   {
     defaultmapfilename = Calibroot;
     defaultmapfilename += "/HCALOUT/tilemap/ohcalgdmlmapfiles102022.root";
   }
  set_default_string_param("MapFileName", defaultmapfilename);
  set_default_string_param("MapHistoName", "ohcal_mephi_map_towerid_");
  
  if (!Calibroot)
  {
    std::cout<<__PRETTY_FUNCTION__ << ": no CALIBRATIONROOT environment variable" << std::endl;
    exit(1);
  }
  set_default_string_param("IronFieldMapPath", std::string(Calibroot) + "/Field/Map/sphenix3dbigmapxyz_steel_rebuild.root" );
  set_default_double_param("IronFieldMapScale", 1.);
  
}
