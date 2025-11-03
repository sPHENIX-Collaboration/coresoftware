#ifndef MICROMEGAS_MicromegasCombinedDataCalibration_H
#define MICROMEGAS_MicromegasCombinedDataCalibration_H

/*!
 * \file MicromegasCombinedDataCalibration.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasMapping.h"

#include <fun4all/SubsysReco.h>

#include <map>
#include <memory>
#include <string>

class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TProfile;

/// micromegas raw data decoder
class MicromegasCombinedDataCalibration : public SubsysReco
{
  public:

  /// constructor
  MicromegasCombinedDataCalibration( const std::string &name = "MicromegasCombinedDataCalibration" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// set min sample for noise estimation
  void set_sample_min( uint16_t value ) { m_sample_min = value; }

  /// set min sample for noise estimation
  void set_sample_max( uint16_t value ) { m_sample_max = value; }

  /// set to true to store evaluation histograms and ntuples
  void set_calibration_file( const std::string& value ) { m_calibration_filename = value; }

  /// do evaluation
  void set_do_evaluation( bool value ) { m_do_evaluation = value; }

  /// evaluation filename
  void set_evaluation_file( const std::string& value ) { m_evaluation_filename = value; }

  private:
  //! raw node
  std::string m_rawhitnodename = "MICROMEGASRAWHIT";

  //! mapping
  MicromegasMapping m_mapping;

  /// min sample for noise estimation
  uint16_t m_sample_min = 0;

  /// max sample for noise estimation
  uint16_t m_sample_max = 100;

  /// calibration output file
  std::string m_calibration_filename = "TPOT_Pedestal_000.root";

  /// map fee id to 2D histograms
  using histogram_map_t = std::map<int, TH2*>;
  histogram_map_t m_histogram_map;

  /// map fee id to Profile histogram
  using profile_map_t = std::map<int, TProfile*>;
  profile_map_t m_profile_map;

  /// do evaluation
  bool m_do_evaluation = false;

  /// evaluation file
  std::string m_evaluation_filename = "MicromegasCombinedDataEvaluation.root";

  // evaluation TFile
  std::unique_ptr<TFile> m_evaluation_file;


};

#endif
