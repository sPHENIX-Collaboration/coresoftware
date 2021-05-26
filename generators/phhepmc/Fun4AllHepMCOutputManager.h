#ifndef PHHEPMC_FUN4ALLHEPMCOUTPUTMANAGER_H
#define PHHEPMC_FUN4ALLHEPMCOUTPUTMANAGER_H

#include <fun4all/Fun4AllOutputManager.h>

#include <fstream>
#include <string>

namespace HepMC
{
  class IO_GenEvent;
}

class PHCompositeNode;

class Fun4AllHepMCOutputManager : public Fun4AllOutputManager
{
 public:
  Fun4AllHepMCOutputManager(const std::string &myname = "HEPMCOUT",
                            const std::string &filename = "hepmcout.txt");

  ~Fun4AllHepMCOutputManager() override;

  int outfileopen(const std::string & /*fname*/) override { return 0; }

  void Print(const std::string &what = "ALL") const override;

  int Write(PHCompositeNode *startNode) override;

  int AddComment(const std::string &text);

  //! embedding ID for the sub-event to be output
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return _embedding_id; }
  //
  //! embedding ID for the sub-event to be output
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { _embedding_id = id; }

 protected:
  std::string outfilename;
  HepMC::IO_GenEvent *ascii_out;
  std::string comment;
  int comment_written;

  // some pointers for use in compression handling
  std::ofstream *filestream;  // holds compressed filestream
  std::ostream *zipstream;    // feed into HepMC

  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int _embedding_id;
};

#endif /* PHHEPMC_FUN4ALLHEPMCOUTPUTMANAGER_H */
