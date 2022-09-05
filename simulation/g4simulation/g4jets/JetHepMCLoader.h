// $Id: $

/*!
 * \file JetHepMCLoader.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4JET_JETHEPMCLOADER_H
#define G4JET_JETHEPMCLOADER_H

#include "Jet.h"

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>

class Fun4AllHistoManager;
class PHCompositeNode;

/*!
 * \brief JetHepMCLoader loads special jet objects encoded in HepMC records to DST Jet nodes. Example use are loading sHijing HIJFRG jets
 *
 * Example use for readback HIJFRAG truth jets from the sHijing HepMC records:
 *
 * \code{.cpp}

    JetHepMCLoader * hepmcjet = new JetHepMCLoader("sHijing_HIJFRG");

    hepmcjet->saveQAPlots();
    hepmcjet->addJet("AntiKt_sHijing_HIJFRG_r02",0,Jet::ANTIKT,0.2,2000000,103);
    hepmcjet->addJet("AntiKt_sHijing_HIJFRG_r04",0,Jet::ANTIKT,0.4,4000000,103);
    hepmcjet->addJet("AntiKt_sHijing_HIJFRG_r06",0,Jet::ANTIKT,0.6,6000000,103);

    se->registerSubsystem(hepmcjet);

 * \endcode
 *
 */
class JetHepMCLoader : public SubsysReco
{
 public:
  //! \param[in] jetInputCategory is the DST PHCompositeNode name that list the output jet maps, e.g. sHijing_HIJFRG for sHijing HIJFRG truth jets
  JetHepMCLoader(const std::string &jetInputCategory);
  ~JetHepMCLoader() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  //! \brief addJet add HepMC entries for a particular type of jets
  //! Example of adding sHijing HIJFRG truth jets with R=0.4:
  //!   \code{.cpp}
  //!     addJet("AntiKt_sHijing_HIJFRG_r04",0,"ANTIKT",0.4,4000000,103);
  //!   \endcode
  //! \param[in] name name of the jet category
  //! \param[in] embeddingID hepmc event's embedding ID
  //! \param[in] algorithm pick one from Jet::ALGO
  //! \param[in] parameter jet parameter, e.g. radius
  //! \param[in] tagPID HepMC entry identifying tag on PID
  //! \param[in] tagStatus HepMC entry identifying tag on status
  void addJet(
      const std::string &name,
      int embeddingID,
      Jet::ALGO algorithm,
      double parameter,
      int tagPID,
      int tagStatus);

  void saveQAPlots(bool b = true) { m_saveQAPlots = b; }

 private:
  int CreateNodes(PHCompositeNode *topNode);
  Fun4AllHistoManager *getHistoManager();

  std::string m_jetInputCategory;

  bool m_saveQAPlots = false;

  struct hepmc_jet_src
  {
    //! name
    std::string m_name;

    //! hepmc event's embedding ID
    int m_embeddingID;

    //! Name of jet algorithm
    std::string m_algorithmName;

    Jet::ALGO m_algorithmID;

    //! jet parameter, e.g. radius
    double m_parameter;

    //! HepMC entry identifying tag on PID
    int m_tagPID;

    //! HepMC entry identifying tag on status
    int m_tagStatus;
  };

  std::vector<hepmc_jet_src> m_jetSrc;
};

#endif /* JETHEPMCLOADER_H_ */
