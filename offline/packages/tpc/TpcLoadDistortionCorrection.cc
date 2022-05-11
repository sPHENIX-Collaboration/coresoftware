/*!
 * \file TpcLoadDistortionCorrection.cc
 * \brief loads distortion correction histogram from file to DistortionCorrectionObject and stores on node tree
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcLoadDistortionCorrection.h"
#include "TpcDistortionCorrectionContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>               

#include <TFile.h>
#include <TH1.h>

namespace
{

  // print histogram
  void print_histogram( TH1* h )
  {

    std::cout << "TpcLoadDistortionCorrection::InitRun - name: " << h->GetName() << std::endl;
    for( const auto& axis:{h->GetXaxis(), h->GetYaxis(), h->GetZaxis() } )
    {
      if( axis ) 
      {
        std::cout
          << "  " << axis->GetName()
          << " bins: " << axis->GetNbins()
          << " min: " << axis->GetXmin()
          << " max: " << axis->GetXmax()
          << std::endl;
      }
    }
    std::cout << std::endl;
  }
  
}  // namespace


//_____________________________________________________________________
TpcLoadDistortionCorrection::TpcLoadDistortionCorrection( const std::string& name ):
  SubsysReco( name)
  {}

//_____________________________________________________________________
int TpcLoadDistortionCorrection::InitRun(PHCompositeNode* topNode)
{

  // look for distortion calibration object
  PHNodeIterator iter(topNode);
 
  /// Get the DST node and check
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "TpcLoadDistortionCorrection::InitRun - DST Node missing, quitting" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
 
  // Get the tracking subnode and create if not found
  auto svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  //create and populate the nodes for each distortion, if present:
  for (int i=0;i<3;i++){

    if( !m_correction_in_use[i] )  continue;

    // get distortion correction object and create if not found
    auto distortion_correction_object = findNode::getClass<TpcDistortionCorrectionContainer>( topNode, m_node_name[i] );
    if( !distortion_correction_object )
      { 
	std::cout << "TpcLoadDistortionCorrection::InitRun - creating TpcDistortionCorrectionContainer in node " << m_node_name[i] << std::endl;
	distortion_correction_object = new TpcDistortionCorrectionContainer;
	auto node = new PHDataNode<TpcDistortionCorrectionContainer>(distortion_correction_object, m_node_name[i]);
	svtxNode->addNode(node);
      }

    std::cout << "TpcLoadDistortionCorrection::InitRun - reading corrections from " << m_correction_filename[i] << std::endl;
    auto distortion_tfile = TFile::Open( m_correction_filename[i].c_str());
    if( !distortion_tfile && m_correction_in_use[i])
      {
	std::cout << "TpcLoadDistortionCorrection::InitRun - cannot open " << m_correction_filename[i] << std::endl;
	exit(1);
      }

    const std::array<const std::string,2> extension = {{ "_negz", "_posz" }};
    for( int i =0; i < 2; ++i )
      {
	distortion_correction_object->m_hDPint[i] = dynamic_cast<TH1*>(distortion_tfile->Get(Form("hIntDistortionP%s", extension[i].c_str()))); assert( distortion_correction_object->m_hDPint[i] );
	distortion_correction_object->m_hDRint[i] = dynamic_cast<TH1*>(distortion_tfile->Get(Form("hIntDistortionR%s", extension[i].c_str()))); assert( distortion_correction_object->m_hDRint[i] );
	distortion_correction_object->m_hDZint[i] = dynamic_cast<TH1*>(distortion_tfile->Get(Form("hIntDistortionZ%s", extension[i].c_str()))); assert( distortion_correction_object->m_hDZint[i] );
      }

      // assign correction object dimension from histograms dimention, assuming all histograms have the same
      distortion_correction_object->dimensions = distortion_correction_object->m_hDPint[0]->GetDimension();
      
      // only dimensions 2 or 3 are supported
      assert( distortion_correction_object->dimensions == 2 || distortion_correction_object->dimensions == 3 );
      
    if( Verbosity() )
    {
      for( const auto& h:{
        distortion_correction_object->m_hDPint[0], distortion_correction_object->m_hDPint[1],
        distortion_correction_object->m_hDRint[0], distortion_correction_object->m_hDRint[1],
        distortion_correction_object->m_hDZint[0], distortion_correction_object->m_hDZint[1] } )
      { print_histogram( h ); }
    }

  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcLoadDistortionCorrection::process_event(PHCompositeNode*)
{ return Fun4AllReturnCodes::EVENT_OK; }
