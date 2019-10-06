// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLPRDFOUTPUTMANAGER_H
#define FUN4ALLRAW_FUN4ALLPRDFOUTPUTMANAGER_H

#include <fun4all/Fun4AllOutputManager.h>

#include <string>

class PHRawOManager;
class PHCompositeNode;

class Fun4AllPrdfOutputManager : public Fun4AllOutputManager
{
 public:
  //! constructor
  Fun4AllPrdfOutputManager(const std::string &myname = "PRDFOUT", const std::string &filename = "data_out.prdf");

  //! destructor
  virtual ~Fun4AllPrdfOutputManager();

  //! PRDF node initialization [class specific method]
  int InitPrdfNode(PHCompositeNode *top_node, const std::string &nodeName = "SIMPRDF");

  //! reinitialize raw output manager to write to new filename. Close old one if any
  int outfileopen(const std::string &fname);

  //! event write method (startNode argument is ignored. prdfNode is always used)
  int Write(PHCompositeNode *startNode);

 private:
  /*! 
    initialize prdf output manager every time 
    the output file name is changed including first event
  */
  int InitPrdfManager();

  //! prdf node
  PHCompositeNode *m_PrdfNode;

  //! output manager
  PHRawOManager *m_PrdfOutManager;
};

#endif
