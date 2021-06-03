#ifndef TPCCALIB_TPCSPACECHARGEMATRIXINVERSION_H
#define TPCCALIB_TPCSPACECHARGEMATRIXINVERSION_H
/**
 * \file TpcSpaceChargeMatrixInversion.h
 * \brief aggregate space charge reconstruction matrices from several jobs and perform the inversion
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */
#include <fun4all/Fun4AllBase.h>

#include <memory>

// forward declaration
class TpcSpaceChargeMatrixContainer;

/**
 * \class TpcSpaceChargeMatrixInversion
 * \brief performs space charge distortion reconstruction using tracks
 */

class TpcSpaceChargeMatrixInversion: public Fun4AllBase
{
  public:

  /// constructor
  TpcSpaceChargeMatrixInversion( const std::string& = "TPCSPACECHARGEMATRIXINVERSION" );

  ///@name modifiers
  //@{

  /// set whether to use only tracks with micromegas or not
  void set_use_micromegas( bool value )
  { m_use_micromegas = value; }

  /// output file
  /**
   * this is the file where space charge correction 3D histograms are stored
   * they are suitable for being read by TpcClusterizer
   */
  void set_outputfile( const std::string& filename );
  
  /// add space charge correction matrix to current. Returns true on success
  bool add( const TpcSpaceChargeMatrixContainer& );

  /// add space charge correction matrix, loaded from file, to current. Returns true on success
  bool add_from_file( const std::string& filename, const std::string& objectname = "TpcSpaceChargeMatrixContainer" );
  
  /// calculate distortions by inverting stored matrices, and save relevant histograms
  void calculate_distortions();

  //@}

  private:

  /// output file
  std::string m_outputfile = "DistortionCorrections.root";

  /// true if only tracks with micromegas must be used
  bool m_use_micromegas = true;

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;

};

#endif
