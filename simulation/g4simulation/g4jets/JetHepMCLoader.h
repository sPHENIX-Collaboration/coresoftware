// $Id: $

/*!
 * \file JetHepMCLoader.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef JETHEPMCLOADER_H_
#define JETHEPMCLOADER_H_

#include <fun4all/SubsysReco.h>
#include <set>
#include <string>

class PHCompositeNode;

/*!
 * \brief JetHepMCLoader loads special jet objects encoded in HepMC records to DST Jet nodes. Example use are loading sHijing HIJFRG jets
 */
class JetHepMCLoader : public SubsysReco
{
 public:

  //! \param[in] jetInputCategory is the DST PHCompositeNode name that list the output jet maps, e.g. sHijing_HIJFRG for sHijing HIJFRG truth jets
  JetHepMCLoader(const std::string &jetInputCategory);
  virtual ~JetHepMCLoader();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  //! \brief addJet add HepMC entries for a particular type of jets
  //! Example of adding sHijing HIJFRG truth jets with R=0.4:
  //!   \code{.cpp}
  //!     addJet("AntiKt_sHijing_HIJFRG_r04",0,"ANTIKT",0.4,4000000,103);
  //!   \endcode
  //! \param[in] name name of the jet category
  //! \param[in] embeddingID hepmc event's embedding ID
  //! \param[in] parameter jet parameter, e.g. radius
  //! \param[in] tagPID HepMC entry identifying tag on PID
  //! \param[in] tagStatus HepMC entry identifying tag on status
  void addJet(
      const std::string &name,
      int embeddingID,
      const std::string &algorithmName,
      double parameter,
      int tagPID,
      int tagStatus);

  void saveQAPlots(bool b) { m_saveQAPlots = b; }

 private:
  std::string m_jetInputCategory;

  bool m_saveQAPlots;

  struct hepmc_jet_src
  {
    //! name
    std::string m_name;

    //! hepmc event's embedding ID
    int m_embeddingID;

    //! Name of jet algorithm
    std::string m_algorithmName;

    //! jet parameter, e.g. radius
    double m_parameter;

    //! HepMC entry identifying tag on PID
    int m_tagPID;

    //! HepMC entry identifying tag on status
    int m_tagStatus;
  };

  int CreateNodes(PHCompositeNode *topNode);

  std::set< hepmc_jet_src> m_jetSrc;

};

#endif /* JETHEPMCLOADER_H_ */
