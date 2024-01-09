#include "CentralityReco.h"

<<<<<<< HEAD
#include "CentralityInfov2.h"

#include <mbd/MbdDefs.h>
#include <mbd/MbdOutV1.h>

#include <fun4all/Fun4AllHistoManager.h>
=======
#include "CentralityInfov1.h"

#include <mbd/MbdOut.h>

#include <ffamodules/CDBInterface.h>

#include <cdbobjects/CDBTTree.h>

>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
<<<<<<< HEAD
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>

#include <phool/phool.h>
#include <phool/recoConsts.h>
=======
#include <phool/getClass.h>
#include <phool/phool.h>

#include <filesystem>
#include <iostream>
#include <string>
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

CentralityReco::CentralityReco(const std::string &name)
  : SubsysReco(name)
{
<<<<<<< HEAD
  _rc = recoConsts::instance();
  _cdb = CDBInterface::instance();
}

CentralityReco::~CentralityReco()
{

}

int CentralityReco::Init(PHCompositeNode * /*unused*/)
{

  _cdb = CDBInterface::instance();

  std::string centdiv_url = _cdb->getUrl("Centrality");

  if (Download_centralityDivisions(centdiv_url)) return Fun4AllReturnCodes::ABORTRUN;

  std::string centscale_url = _cdb->getUrl("CentralityScale");

  if (Download_centralityScale(centscale_url)) return Fun4AllReturnCodes::ABORTRUN;

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  CreateNodes(topNode);
  return 0;
}

int CentralityReco::Download_centralityScale(const std::string& dbfile)
{

  _centrality_scale = 1.00;

  TString dbase_file = dbfile;

  if (dbase_file.EndsWith(".root"))
    {
      CDBTTree *cdbttree = new CDBTTree(dbase_file.Data());
      cdbttree->LoadCalibrations();
      _centrality_scale = cdbttree->GetDoubleValue(0,"centralityscale");
      if (Verbosity()) std::cout << "centscale = "<<_centrality_scale << std::endl;
      delete cdbttree;
    }
  else 
    {
      std::cout << PHWHERE <<", ERROR, unknown file type, " << dbfile <<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
=======
}

int CentralityReco::InitRun(PHCompositeNode *topNode)
{
  CDBInterface *_cdb = CDBInterface::instance();

  std::string centdiv_url = _cdb->getUrl("Centrality");

  if (Download_centralityDivisions(centdiv_url))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::string centscale_url = _cdb->getUrl("CentralityScale");

  if (Download_centralityScale(centscale_url))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CreateNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::Download_centralityScale(const std::string &dbfile)
{
  _centrality_scale = 1.00;

  std::filesystem::path dbase_file = dbfile;

  if (dbase_file.extension() == ".root")
  {
    CDBTTree *cdbttree = new CDBTTree(dbase_file);
    cdbttree->LoadCalibrations();
    _centrality_scale = cdbttree->GetDoubleValue(0, "centralityscale");
    if (Verbosity())
    {
      std::cout << "centscale = " << _centrality_scale << std::endl;
    }
    delete cdbttree;
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbfile << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  return Fun4AllReturnCodes::EVENT_OK;
}

<<<<<<< HEAD
int CentralityReco::Download_centralityDivisions(const std::string& dbfile)
{

  for (int idiv = 0; idiv < 20; idiv++)
    _centrality_map[idiv] = 0;

  TString dbase_file = dbfile;

  if (dbase_file.EndsWith(".root"))
    {
      CDBTTree *cdbttree = new CDBTTree(dbase_file.Data());
      cdbttree->LoadCalibrations();
      for (int idiv = 0; idiv < NDIVS;idiv++)
	{
	  _centrality_map[idiv] = cdbttree->GetFloatValue(idiv,"centralitydiv");
	  if (Verbosity()) std::cout << "centdiv "<<idiv<<" : "<<_centrality_map[idiv]<<std::endl;
	}
      delete cdbttree;
    }
  else 
    {
      std::cout << PHWHERE <<", ERROR, unknown file type, " << dbfile <<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
=======
int CentralityReco::Download_centralityDivisions(const std::string &dbfile)
{
  _centrality_map.fill(0);

  std::filesystem::path dbase_file = dbfile;

  if (dbase_file.extension() == ".root")
  {
    CDBTTree *cdbttree = new CDBTTree(dbase_file);
    cdbttree->LoadCalibrations();
    for (int idiv = 0; idiv < NDIVS; idiv++)
    {
      _centrality_map[idiv] = cdbttree->GetFloatValue(idiv, "centralitydiv");
      if (Verbosity())
      {
        std::cout << "centdiv " << idiv << " : " << _centrality_map[idiv] << std::endl;
      }
    }
    delete cdbttree;
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbfile << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  return Fun4AllReturnCodes::EVENT_OK;
}

<<<<<<< HEAD

int CentralityReco::ResetEvent(PHCompositeNode *)
=======
int CentralityReco::ResetEvent(PHCompositeNode * /*unused*/)
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  _mbd_charge_sum = 0.;
  _mbd_charge_sum_n = 0.;
  _mbd_charge_sum_s = 0.;

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::FillVars()
{
<<<<<<< HEAD

  if (Verbosity()>1)
=======
  if (Verbosity() > 1)
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  _mbd_charge_sum_s = _mbd_out->get_q(0);

  _mbd_charge_sum_n = _mbd_out->get_q(1);

<<<<<<< HEAD
  _mbd_charge_sum = (_mbd_charge_sum_n + _mbd_charge_sum_s)*_centrality_scale;


  if (Verbosity())
  {
    std::cout << "  MBD sum = " <<_mbd_charge_sum<<std::endl;
=======
  _mbd_charge_sum = (_mbd_charge_sum_n + _mbd_charge_sum_s) * _centrality_scale;

  if (Verbosity())
  {
    std::cout << "  MBD sum = " << _mbd_charge_sum << std::endl;
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
    std::cout << "      North: " << _mbd_charge_sum_n << std::endl;
    std::cout << "      South: " << _mbd_charge_sum_s << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

<<<<<<< HEAD

int CentralityReco::FillCentralityInfo()
{
  // Fill is minbias
  

  if (Verbosity()>1)
=======
int CentralityReco::FillCentralityInfo()
{
  // Fill is minbias

  if (Verbosity() > 1)
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

<<<<<<< HEAD

  float value = -999.99;
=======
  float value = std::numeric_limits<float>::quiet_NaN();
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  for (int i = 0; i < NDIVS; i++)
  {
    if (_centrality_map[i] < _mbd_charge_sum)
    {
<<<<<<< HEAD

=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
      value = 0.05 * i;
      break;
    }
  }

  _central->set_centile(CentralityInfo::PROP::mbd_NS, value);

  return Fun4AllReturnCodes::EVENT_OK;
<<<<<<< HEAD

=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
}

int CentralityReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << "------------CentralityReco-------------" << std::endl;
  }

<<<<<<< HEAD

=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  // Get Nodes from the Tree
  if (GetNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Fill Arrays
  if (FillVars())
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (FillCentralityInfo())
<<<<<<< HEAD
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  if (Verbosity())
    {
      _central->identify();
    }
=======
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (Verbosity())
  {
    _central->identify();
  }
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

<<<<<<< HEAD
  _central = findNode::getClass<CentralityInfov2>(topNode, "CentralityInfo");
  
  if (!_central)
    {
      std::cout << "no centrality node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    
  
  _mbd_out = findNode::getClass<MbdOutV1>(topNode, "MbdOut");
  
  if (!_mbd_out)
    {
      std::cout << "no MBD out node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
=======
  _central = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");

  if (!_central)
  {
    std::cout << "no centrality node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _mbd_out = findNode::getClass<MbdOut>(topNode, "MbdOut");

  if (!_mbd_out)
  {
    std::cout << "no MBD out node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  return Fun4AllReturnCodes::EVENT_OK;
}

void CentralityReco::CreateNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
  }

  PHNodeIterator dstIter(dstNode);

  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!detNode)
  {
    std::cout << PHWHERE << "Detector Node missing, making one" << std::endl;
    detNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(detNode);
  }

<<<<<<< HEAD
  CentralityInfov2 *central = new CentralityInfov2();
  
=======
  CentralityInfo *central = new CentralityInfov1();

>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  PHIODataNode<PHObject> *centralityNode = new PHIODataNode<PHObject>(central, "CentralityInfo", "PHObject");
  detNode->addNode(centralityNode);

  return;
}
<<<<<<< HEAD

int CentralityReco::End(PHCompositeNode * /* topNode*/)
{

  return 0;
}
=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
