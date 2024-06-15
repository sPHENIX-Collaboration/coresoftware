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
    DistortionType_Fluctuation = 2
  };

  //! correction filename
  void set_correction_filename(DistortionType i, const std::string& value)
  {
    if (i < 0 || i >= 3) return;
    m_correction_filename[i] = value;
    m_correction_in_use[i] = true;
  }

  //! set the phi histogram to be interpreted as radians.
  void set_read_phi_as_radians(bool flag)
  {
    m_phi_hist_in_radians = flag;
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
  std::string m_correction_filename[3] = {"", "", ""};

  //! flag to indicate correction in use
  bool m_correction_in_use[3] = {false, false, false};

  //! set the phi histogram to be interpreted as radians rather than mm
  bool m_phi_hist_in_radians = true;

  //! distortion object node name
  std::string m_node_name[3] = {"TpcDistortionCorrectionContainerStatic", "TpcDistortionCorrectionContainerAverage", "TpcDistortionCorrectionContainerFluctuation"};
};

#endif
