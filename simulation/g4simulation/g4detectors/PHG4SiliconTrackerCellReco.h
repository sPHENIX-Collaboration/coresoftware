#ifndef PHG4SILICONTRACKERCELLRECO_H
#define PHG4SILICONTRACKERCELLRECO_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#ifndef __CINT__
#include <gsl/gsl_vector.h>
#endif

#include <string>

class PHCompositeNode;
class PHG4Cell;

class PHG4SiliconTrackerCellReco : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4SiliconTrackerCellReco(const std::string &name = "SILICON_TRACKER");

  virtual ~PHG4SiliconTrackerCellReco();
  //! module initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! set default parameter values
  void SetDefaultParameters();

  void Detector(const std::string &d) { m_Detector = d; }
  void checkenergy(const int i = 1) { m_ChkEnergyConservationFlag = i; }

  double circle_rectangle_intersection(double x1, double y1,  double x2,  double y2,  double mx, double my,  double r) const;
  double sA(double r, double x, double y) const;
  
 protected:
  int CheckEnergy(PHCompositeNode *topNode);
  std::string m_Detector;
  std::string m_HitNodeName;
  std::string m_CellNodeName;
  std::string m_GeoNodeName;

  int m_ChkEnergyConservationFlag;
  std::map<std::string, PHG4Cell *> m_CellList;  // This map holds the hit cells

  double m_Tmin;
  double m_Tmax;

#ifndef __CINT__
  gsl_vector *m_LocalOutVec;
  gsl_vector *m_PathVec;
  gsl_vector *m_SegmentVec;
#endif
};

#endif
