#ifndef TPC_TPCLOADDISTORTIONCORRECTION_H
#define TPC_TPCLOADDISTORTIONCORRECTION_H

/*!
 * \file TpcLoadDistortionCorrection.h
 * \brief loads distortion correction histogram from file to DistortionCorrectionObject and stores on node tree
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <trackbase/TrkrDefs.h>

#include <array>

class TH3;

class TpcLoadDistortionCorrection : public SubsysReco
{
 public:
  //! constructor
  TpcLoadDistortionCorrection(const std::string& = "TpcLoadDistortionCorrection");

  //! global initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  // convenient enumeration for distortion type
  enum DistortionType:int
  {
    DistortionType_Static = 0,
    DistortionType_Average = 1,
    DistortionType_Fluctuation = 2,
    DistortionType_ModuleEdge = 3
  };

  static const int nDistortionTypes = 4;

  //! correction filename
  void set_correction_filename(DistortionType i, const std::string& value)
  {
    if (i < 0 || i >= nDistortionTypes) return;
    m_correction_filename[i] = value;
    m_correction_in_use[i] = true;
  }

  //! set the phi histogram to be interpreted as radians.
  void set_read_phi_as_radians(bool flag)
  {
    m_phi_hist_in_radians[0] = flag;
  }
 void set_read_phi_as_radians(int i, bool flag)
  {
    m_phi_hist_in_radians[i] = flag;
  }
  //! set the histogram to interpolate between hist value and zero, depending on z position. (has no effect if m_dimensions is 3)
  void set_interpolate_2D_to_zero(bool flag)
  {
    m_interpolate_z[0] = flag;
  }
  void set_interpolate_2D_to_zero(int i, bool flag)
  {
    m_interpolate_z[i] = flag;
  } 

  //! node name
  void set_node_name(const std::string& value)
  {
    m_node_name[0] = value;
  }
  void set_node_name(int i, const std::string& value)
  {
    m_node_name[i] = value;
  }

 private:
  //! correction filename
  std::string m_correction_filename[nDistortionTypes] = {"", "", "",""};

  //! flag to indicate correction in use
  bool m_correction_in_use[nDistortionTypes] = {false, false, false,false};

  //! set the phi histogram to be interpreted as radians rather than mm
  bool m_phi_hist_in_radians[nDistortionTypes] = {true,true,true,true};
  bool m_interpolate_z[nDistortionTypes] = {true,true,true,true};

  //! distortion object node name
  std::string m_node_name[nDistortionTypes] = {"TpcDistortionCorrectionContainerStatic", "TpcDistortionCorrectionContainerAverage", "TpcDistortionCorrectionContainerFluctuation","TpcDistortionCorrectionContainerModuleEdge"};
};

#endif
