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
  TpcLoadDistortionCorrection(  const std::string& = "TpcLoadDistortionCorrection" );

  //! global initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  //! distortion filename
  void set_distortion_filename( const std::string& value )
  { m_distortion_filename = value; }

  //! node name
  void set_node_name( const std::string& value )
  { m_node_name = value; }
  
  private:

  //! distortion filename
  std::string m_distortion_filename;

  //! distortion object node name
  std::string m_node_name = "TpcDistortionCorrectionContainer";

};

#endif
