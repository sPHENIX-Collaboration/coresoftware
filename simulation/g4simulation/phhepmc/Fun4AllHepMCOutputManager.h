#ifndef FUN4ALLHEPMCOUTPUTMANAGER_H__
#define FUN4ALLHEPMCOUTPUTMANAGER_H__


#include <fun4all/Fun4AllOutputManager.h>
#include <string>

#include <fstream>
#include <ostream>

  namespace HepMC {
    class IO_GenEvent;
  };


class PHCompositeNode;

class Fun4AllHepMCOutputManager: public Fun4AllOutputManager
{
 public:

  Fun4AllHepMCOutputManager(const std::string &myname = "HEPMCOUT", 
			    const std::string &filename = "hepmcout.txt");

  virtual ~Fun4AllHepMCOutputManager();

  int outfileopen(const std::string& /*fname*/) {return 0;}

  void Print(const std::string &what = "ALL") const;

  int Write(PHCompositeNode *startNode);

  int AddComment(const std::string &text);

 protected:
  std::string outfilename;
  HepMC::IO_GenEvent *ascii_out;
  std::string comment;
  int comment_written;

  // some pointers for use in compression handling
  std::ofstream *filestream; // holds compressed filestream
  std::ostream *zipstream;   // feed into HepMC
};

#endif /* FUN4ALLHEPMCOUTPUTMANAGER_H__ */
