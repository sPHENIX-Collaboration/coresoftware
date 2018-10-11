#include "Fun4AllPrdfOutputManager.h"
#include <phool/recoConsts.h>

#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRawOManager.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

//______________________________________________________
Fun4AllPrdfOutputManager::Fun4AllPrdfOutputManager( const string &myname, const string &fname):
  Fun4AllOutputManager( myname ),
  prdfNode( 0 ),
  prdfOut( 0 )
{
  outfilename = fname;
  return ;
}

//______________________________________________________
Fun4AllPrdfOutputManager::~Fun4AllPrdfOutputManager( void )
{
  if( prdfOut ) delete prdfOut;
  return;
}
  
//______________________________________________________
int Fun4AllPrdfOutputManager::InitPrdfNode( PHCompositeNode* top_node, const string &nodeName )
{
  PHNodeIterator topIter(top_node);
  prdfNode = dynamic_cast<PHCompositeNode*>(topIter.findFirst("PHCompositeNode", nodeName.c_str()));
  if ( prdfNode )
    {
      // the prdfNode already exists (Pisa Input Mgr creates one also)
      return 0;
    }

  // check name wrt default
  if ( nodeName != "SIMPRDF" )
    cout << "Fun4AllPrdfOutputManager::InitPrdfNode - WARNING: nodeName is \"" << nodeName << "\". most systems expect \"SIMPRDF\" and this is most likely not going to work" << endl;

  // create node
  prdfNode = new PHCompositeNode( nodeName.c_str() );
  top_node->addNode(prdfNode);
  return 0;

}
  
//______________________________________________________
int Fun4AllPrdfOutputManager::outfileopen(const string &fname)
{
  if( prdfOut ) {
    if( Verbosity() ) cout << "Fun4AllPrdfOutputManager::outfileopen - closing file \"" << outfilename << "\"" << endl;
    delete prdfOut;
    prdfOut = 0;
  }

  outfilename = fname;
  if( Verbosity() ) cout << "Fun4AllPrdfOutputManager::outfileopen - writing to file \"" << outfilename << "\"" << endl;

  return 0;
  
}

//______________________________________________________
int Fun4AllPrdfOutputManager::Write(PHCompositeNode* /*node*/)
{
  
  // check prdfNode
  if( !prdfNode ) {
    cout << "Fun4AllPrdfOutputManager::Write - prdfNode not initialized" << endl;
    return -1;
  }
  
  // check prdfOut
  if( !prdfOut ) InitPrdfManager();
  if( !prdfOut ) {
    cout << "Fun4AllPrdfOutputManager::Write - prdf manager not initialized" << endl;
    return -1;
  }

  // write prdfNode to prdfManager
  bool prdf_status = prdfOut->write(prdfNode);
  return prdf_status ? 0:-1;

}

//______________________________________________________
int Fun4AllPrdfOutputManager::InitPrdfManager( void )
{

  if( prdfOut ) return -1;
  
  // retrieve run number from recoConsts
  recoConsts *rc = recoConsts::instance();
  
  // retrieve run number
  int run_number = -1;
  if( rc->FlagExist( "RUNNUMBER") )
    { 
      run_number = rc->get_IntFlag("RUNNUMBER");
    }
  // buffer length (taken from offline/framework/simReco/PrdfReco.C)
  static const int buffer_length( 8*1024*1024/4 );
  
  // create output manager
  prdfOut = new PHRawOManager( outfilename.c_str() , run_number , buffer_length);
  return 0;

}
