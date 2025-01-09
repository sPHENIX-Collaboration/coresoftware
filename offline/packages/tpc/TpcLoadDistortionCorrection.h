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

  static constexpr int nDistortionTypes = 4;

  //! correction filename
  void set_correction_filename(DistortionType i, const std::string& value)
  {
    if (i < 0 || i >= nDistortionTypes) return;
    m_correction_filename[i] = value;
    m_correction_in_use[i] = true;
  }

  //! set the scale factor to be applied to the correction
  void set_scale_factor(DistortionType i, float value)
  {
    m_use_scalefactor[i] = true;
    m_scalefactor[i] = value;
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
  // false:  correction is corr(x,y) (regardless of z)
  // true:  correction is corr(x,y)*(1-(z/zspan)) so that corr(x,y) at readout is zero.
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
  std::array<std::string,nDistortionTypes> m_correction_filename = {};

  //! flag to indicate correction in use
  std::array<bool,nDistortionTypes> m_correction_in_use = {};

  //! flag and scalefactor to apply to correction
  std::array<bool,nDistortionTypes> m_use_scalefactor = {};

  //! scale factors
  std::array<float,nDistortionTypes> m_scalefactor = {1.0,1.0,1.0,1.0};

  //! set the phi histogram to be interpreted as radians rather than mm
  std::array<bool,nDistortionTypes> m_phi_hist_in_radians = {true,true,true,true};

  //! z interpolation
  std::array<bool,nDistortionTypes> m_interpolate_z = {true,true,true,true};

  //! distortion object node name
  std::array<std::string,nDistortionTypes> m_node_name = {"TpcDistortionCorrectionContainerStatic", "TpcDistortionCorrectionContainerAverage", "TpcDistortionCorrectionContainerFluctuation","TpcDistortionCorrectionContainerModuleEdge"};
};

#endif
