#ifndef FUN4ALL_FUN4ALLPRDFOUTPUTMANAGER_H
#define FUN4ALL_FUN4ALLPRDFOUTPUTMANAGER_H

#include "Fun4AllOutputManager.h"

#include <string>

class PHRawOManager;
class PHCompositeNode;

class Fun4AllPrdfOutputManager: public Fun4AllOutputManager
{
 public:

  //! constructor
  Fun4AllPrdfOutputManager(const std::string &myname = "PRDFOUT" , const std::string &filename = "data_out.prdf");

  //! destructor
  virtual ~Fun4AllPrdfOutputManager( void );
  
  //! PRDF node initialization [class specific method]
  int InitPrdfNode( PHCompositeNode* top_node, const std::string &nodeName="SIMPRDF" );
  
  //! reinitialize raw output manager to write to new filename. Close old one if any
  int outfileopen(const std::string &fname);
  
  //! event write method (startNode argument is ignored. prdfNode is always used)
  int Write(PHCompositeNode *startNode);

 protected:

  /*! 
    initialize prdf output manager every time 
    the output file name is changed including first event
  */
  int InitPrdfManager( void );

  //! prdf node
  PHCompositeNode *prdfNode;
  
  //! output manager
  PHRawOManager *prdfOut;
};

#endif
  
