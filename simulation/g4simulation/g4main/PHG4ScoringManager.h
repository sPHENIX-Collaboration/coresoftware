// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4ScoringManager.h
 * \brief the connection between Fun4All to G4ScoringManager
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4MAIN_PHG4SCORINGMANAGER_H
#define G4MAIN_PHG4SCORINGMANAGER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <utility>

class Fun4AllHistoManager;
class PHCompositeNode;

/*!
 * \brief PHG4ScoringManager is the connection between Fun4All to G4ScoringManager
 *  Track primitive score like flux or energy deposition integrated over events and save to histograms
 *  More on G4ScoringManager see
 *  Manual http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/BackupVersions/V9.6/html/ch04s08.html
 *  And talk http://geant4.slac.stanford.edu/JLAB2012/Scoring1.pdf
 */
class PHG4ScoringManager : public SubsysReco
{
 public:
  PHG4ScoringManager();

  ~PHG4ScoringManager() override {}

  //! full initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! end of run method
  int End(PHCompositeNode *) override;

  //! Output result to a ROOT file with this name
  void setOutputFileName(const std::string &outputfilename) { m_outputFileName = outputfilename; };

  //! \brief Run this Geant4 command after initialization of G4ScoringManager in the InitRun() stage
  //! You can call this command multiple times to stage multiple commands to run
  /*!
   *  Example
      \code{.cpp}

      PHG4ScoringManager * g4score = new PHG4ScoringManager();
      g4score->G4Command("/score/create/cylinderMesh cyl_score");
      g4score->G4Command("/score/mesh/cylinderSize 100. 400. cm");
      g4score->G4Command("/score/mesh/nBin 100 800 8");
      g4score->G4Command("/score/quantity/eDep cyl_edep");
      g4score->G4Command("/score/quantity/cellFlux cyl_flux");
      g4score->G4Command("/score/close");

      \endcode
   */
  void G4Command(const std::string &cmd);

  void setVertexHistRange(double min, double max) {m_vertexHistRange.first = min; m_vertexHistRange.second = max;}
  void setVertexAcceptanceRange(double min, double max) {m_vertexAcceptanceRange.first = min; m_vertexAcceptanceRange.second = max;}

 private:
  Fun4AllHistoManager *getHistoManager();
  void makeScoringHistograms();

  std::vector<std::string> m_commands;

  std::string m_outputFileName;

  std::pair<double, double> m_vertexHistRange {-5e2 ,5e2 };
  std::pair<double, double> m_vertexAcceptanceRange {-5e2 ,5e2 };
};

#endif /* PHG4SCORINGMANAGER_H_ */
