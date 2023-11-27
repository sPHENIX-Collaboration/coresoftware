#include "MbdCalib.h"
#include <MbdGeomV1.h>

// Database Includes
#include <ffamodules/CDBInterface.h>

#include <cdbobjects/CDBTTree.h>

#include <phool/phool.h>

#include <cmath>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>

MbdCalib::MbdCalib()
{
  Reset();
  _rc = recoConsts::instance();
  _mbdgeom = std::make_unique<MbdGeomV1>();

  if ( _rc->FlagExist("MBD_TEMPLATEFIT") )
  {
    do_templatefit = _rc->get_IntFlag("MBD_TEMPLATEFIT");
  }
  else
  {
    do_templatefit = 0;
  }

}

int MbdCalib::Download_All()
{
  _status = 0;

  if (Verbosity() > 0)
  {
    std::cout << "MBD CDB " << _rc->get_StringFlag("CDB_GLOBALTAG") << "\t" << _rc->get_uint64Flag("TIMESTAMP") << std::endl;
  }
  _cdb = CDBInterface::instance();
  // if rc flag MBD_CALDIR does not exist, we create it and set it to an empty string
  if (!_rc->FlagExist("MBD_CALDIR"))
  {
    std::string sampmax_url = _cdb->getUrl("MBD_SAMPMAX");
    // std::cout << "sampmax_url " << sampmax_url << std::endl;
    Download_SampMax(sampmax_url);

    std::string qfit_url = _cdb->getUrl("MBD_QFIT");
    // std::cout << "qfit_url " << qfit_url << std::endl;
    Download_Gains(qfit_url);

    std::string tq_t0_url = _cdb->getUrl("MBD_TQ_T0");
    // std::cout << "tq_t0_url " << tq_t0_url << std::endl;
    Download_TQT0(tq_t0_url);

    if ( do_templatefit )
    {
      std::string shape_url = _cdb->getUrl("MBD_SHAPES");
      Download_Shapes(shape_url);
    }
  }
  else
  {
    std::string bbc_caldir = _rc->get_StringFlag("MBD_CALDIR");
    std::cout << "Reading MBD Calibrations from " << bbc_caldir << std::endl;

    std::string sampmax_file = bbc_caldir + "/bbc_sampmax.calib";
    Download_SampMax(sampmax_file);
    std::string qfit_file = bbc_caldir + "/bbc_qfit.calib";
    Download_Gains(qfit_file);
    std::string tq_t0_file = bbc_caldir + "/bbc_tq_t0.calib";
    Download_TQT0(tq_t0_file);

    if ( do_templatefit )
    {
      std::string shape_file = bbc_caldir + "/bbc_shape.calib";
      Download_Shapes(shape_file);
    }
  }

  return _status;
}

int MbdCalib::Download_Gains(const std::string& dbase_location)
{
  // Reset All Values
  _qfit_integ.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_mpv.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_sigma.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_integerr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_mpverr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_sigmaerr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_chi2ndf.fill(std::numeric_limits<float>::quiet_NaN());

  if (Verbosity() > 0)
  {
    std::cout << "Opening " << dbase_location << std::endl;
  }
  std::filesystem::path dbase_file = dbase_location;
  if (dbase_file.extension() == ".root")  // read from database
  {
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    cdbttree->LoadCalibrations();

    for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
    {
      _qfit_integ[ipmt] = cdbttree->GetFloatValue(ipmt, "qfit_integ");
      _qfit_mpv[ipmt] = cdbttree->GetFloatValue(ipmt, "qfit_mpv");
      _qfit_sigma[ipmt] = cdbttree->GetFloatValue(ipmt, "qfit_sigma");
      _qfit_integerr[ipmt] = cdbttree->GetFloatValue(ipmt, "qfit_integerr");
      _qfit_mpverr[ipmt] = cdbttree->GetFloatValue(ipmt, "qfit_mpverr");
      _qfit_sigmaerr[ipmt] = cdbttree->GetFloatValue(ipmt, "qfit_sigmaerr");
      _qfit_chi2ndf[ipmt] = cdbttree->GetFloatValue(ipmt, "qfit_chi2ndf");
      if (Verbosity() > 0)
      {
        if (ipmt < 5)
        {
          std::cout << ipmt << "\t" << _qfit_mpv[ipmt] << std::endl;
        }
      }
    }
    delete cdbttree;
  }
  else if (dbase_file.extension() == ".calib")  // read from text file
  {
    std::ifstream infile(dbase_location);
    if (!infile.is_open())
    {
      std::cout << PHWHERE << "unable to open " << dbase_location << std::endl;
      _status = -3;
      return _status;
    }

    int pmt = -1;
    while (infile >> pmt)
    {
      infile >> _qfit_integ[pmt] >> _qfit_mpv[pmt] >> _qfit_sigma[pmt] >> _qfit_integerr[pmt] >> _qfit_mpverr[pmt] >> _qfit_sigmaerr[pmt] >> _qfit_chi2ndf[pmt];
      if (Verbosity() > 0)
      {
        if (pmt < 5 || pmt >= MbdDefs::MBD_N_PMT - 5)
        {
          std::cout << pmt << "\t" << _qfit_integ[pmt] << "\t" << _qfit_mpv[pmt] << "\t" << _qfit_sigma[pmt]
                    << "\t" << _qfit_integerr[pmt] << "\t" << _qfit_mpverr[pmt] << "\t" << _qfit_sigmaerr[pmt]
                    << "\t" << _qfit_chi2ndf[pmt] << std::endl;
        }
      }
    }
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;
  }

  return 1;
}

int MbdCalib::Download_TQT0(const std::string& dbase_location)
{
  // Reset All Values
  _tqfit_t0mean.fill(std::numeric_limits<float>::quiet_NaN());
  _tqfit_t0meanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _tqfit_t0sigma.fill(std::numeric_limits<float>::quiet_NaN());
  _tqfit_t0sigmaerr.fill(std::numeric_limits<float>::quiet_NaN());

  if (Verbosity() > 0)
  {
    std::cout << "Opening " << dbase_location << std::endl;
  }
  std::filesystem::path dbase_file = dbase_location;
  if (dbase_file.extension() == ".root")  // read from database
  {
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    cdbttree->LoadCalibrations();

    for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
    {
      _tqfit_t0mean[ipmt] = cdbttree->GetFloatValue(ipmt, "tqfit_t0mean");
      _tqfit_t0meanerr[ipmt] = cdbttree->GetFloatValue(ipmt, "tqfit_t0meanerr");
      _tqfit_t0sigma[ipmt] = cdbttree->GetFloatValue(ipmt, "tqfit_t0sigma");
      _tqfit_t0sigmaerr[ipmt] = cdbttree->GetFloatValue(ipmt, "tqfit_t0sigmaerr");
      if (Verbosity() > 0)
      {
        if (ipmt < 5 || ipmt >= MbdDefs::MBD_N_PMT - 5)
        {
          std::cout << ipmt << "\t" << _tqfit_t0mean[ipmt] << std::endl;
        }
      }
    }
    delete cdbttree;
  }
  else if (dbase_file.extension() == ".calib")  // read from text file
  {
    std::ifstream infile(dbase_location);
    if (!infile.is_open())
    {
      std::cout << PHWHERE << "unable to open " << dbase_location << std::endl;
      _status = -3;
      return _status;
    }

    int pmt = -1;
    while (infile >> pmt)
    {
      infile >> _tqfit_t0mean[pmt] >> _tqfit_t0meanerr[pmt] >> _tqfit_t0sigma[pmt] >> _tqfit_t0sigmaerr[pmt];

      if (Verbosity() > 0)
      {
        if (pmt < 5 || pmt >= MbdDefs::MBD_N_PMT - 5)
        {
          std::cout << pmt << "\t" << _tqfit_t0mean[pmt] << "\t" << _tqfit_t0meanerr[pmt]
                    << "\t" << _tqfit_t0sigma[pmt] << "\t" << _tqfit_t0sigmaerr[pmt] << std::endl;
        }
      }
    }
    infile.close();
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;
  }

  return 1;
}

int MbdCalib::Download_SampMax(const std::string& dbase_location)
{
  // Reset All Values
  _sampmax.fill(-1);

  std::filesystem::path dbase_file = dbase_location;
  if (dbase_file.extension() == ".root")  // read from database
  {
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    cdbttree->LoadCalibrations();

    for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
    {
      _sampmax[ifeech] = cdbttree->GetIntValue(ifeech, "sampmax");
      if (Verbosity() > 0)
      {
        if (ifeech < 5 || ifeech >= MbdDefs::MBD_N_FEECH - 5)
        {
          std::cout << ifeech << "\t" << _sampmax[ifeech] << std::endl;
        }
      }
    }
    delete cdbttree;
  }
  else if (dbase_file.extension() == ".calib")  // read from text file
  {
    std::ifstream infile(dbase_location);
    if (!infile.is_open())
    {
      std::cout << PHWHERE << "unable to open " << dbase_location << std::endl;
      _status = -3;
      return _status;
    }

    int feech = -1;
    while (infile >> feech)
    {
      infile >> _sampmax[feech];
      if (Verbosity() > 0)
      {
        if (feech < 5 || feech >= MbdDefs::MBD_N_FEECH - 5)
        {
          std::cout << feech << "\t" << _sampmax[feech] << std::endl;
        }
      }
    }
    infile.close();
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;  // file not found
  }

  return 1;
}

int MbdCalib::Download_Shapes(const std::string& dbase_location)
{
  Verbosity(100);
  if ( Verbosity() ) std::cout << "In MbdCalib::Download_Shapes" << std::endl;
  // Reset All Values
  for(auto& shape : _shape_y) {
    shape.clear();
  }
  for(auto& sherr : _sherr_yerr) {
    sherr.clear();
  }

  std::filesystem::path dbase_file = dbase_location;
  if (dbase_file.extension() == ".root")  // read from database
  {
    if ( Verbosity() ) std::cout << "Reading from CDB " << dbase_location << std::endl;
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    cdbttree->LoadCalibrations();

    for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
    {
      if ( _mbdgeom->get_type(ifeech) == 0 ) continue;  // skip t-channels

      _shape_npts[ifeech] = cdbttree->GetIntValue(ifeech, "shape_npts");
      _shape_minrange[ifeech] = cdbttree->GetDoubleValue(ifeech, "shape_min");
      _shape_maxrange[ifeech] = cdbttree->GetDoubleValue(ifeech, "shape_max");

      _sherr_npts[ifeech] = cdbttree->GetIntValue(ifeech, "sherr_npts");
      _sherr_minrange[ifeech] = cdbttree->GetDoubleValue(ifeech, "sherr_min");
      _sherr_maxrange[ifeech] = cdbttree->GetDoubleValue(ifeech, "sherr_max");

      for (int ipt=0; ipt<_shape_npts[ifeech]; ipt++)
      {
        int chtemp = 1000*ipt + ifeech;

        Double_t val = cdbttree->GetDoubleValue(chtemp, "shape_val");
        _shape_y[ifeech].push_back( val );

        val = cdbttree->GetDoubleValue(chtemp, "sherr_val");
        _sherr_yerr[ifeech].push_back( val );
      }

      if (Verbosity() > 0)
      {
        if (ifeech < 5 || ifeech >= MbdDefs::MBD_N_FEECH - 5)
        {
          std::cout << ifeech << "\t" << _shape_y[ifeech][0] << std::endl;
        }
      }
    }
    delete cdbttree;
  }
  else if (dbase_file.extension() == ".calib")  // read from text file
  {
    if ( Verbosity() ) std::cout << "Reading from " << dbase_location << std::endl;
    std::ifstream infile(dbase_location);
    if (!infile.is_open())
    {
      std::cout << PHWHERE << "unable to open " << dbase_location << std::endl;
      _status = -3;
      return _status;
    }

    int temp_feech = -1;
    int temp_npoints = -1;
    float temp_begintime = -1;
    float temp_endtime = -1;
    while ( infile >> temp_feech >> temp_npoints >> temp_begintime >> temp_endtime )
    {
      if ( Verbosity() ) std::cout << "shape " << temp_feech << "\t" <<  temp_npoints << "\t" <<  temp_begintime << "\t" <<  temp_endtime << std::endl;
      if ( temp_feech<0 || temp_feech>255 )
      {
        std::cout << "ERROR, invalid FEECH " << temp_feech << " in MBD waveforms calibration" << std::endl;
        _status = -2;
        return _status;
      }

      _shape_npts[temp_feech] = temp_npoints;
      _shape_minrange[temp_feech] = temp_begintime;
      _shape_maxrange[temp_feech] = temp_endtime;

      Double_t temp_val{0.};
      for (int isamp=0; isamp<temp_npoints; isamp++)
      {
        infile >> temp_val;
        _shape_y[temp_feech].push_back( temp_val );
        if ( Verbosity() && (temp_feech==8 || temp_feech==255) )
        {
          std::cout << _shape_y[temp_feech][isamp] << " ";
          if ( isamp%10==9 ) std::cout << std::endl;
        }
      }
      if ( Verbosity() ) std::cout << std::endl;
    }

    infile.close();

    // Now read in the sherr file
    std::string sherr_dbase_location = std::regex_replace( dbase_location, std::regex("bbc_shape.calib"), "bbc_sherr.calib" );
    if ( Verbosity() ) std::cout << "Reading from " << sherr_dbase_location << std::endl;
    infile.open( sherr_dbase_location );
    if (!infile.is_open())
    {
      std::cout << PHWHERE << "unable to open " << sherr_dbase_location << std::endl;
      _status = -3;
      return _status;
    }

    temp_feech = -1;
    temp_npoints = -1;
    temp_begintime = -1;
    temp_endtime = -1;
    while ( infile >> temp_feech >> temp_npoints >> temp_begintime >> temp_endtime )
    {
      if ( Verbosity() ) std::cout << "sheer " << temp_feech << "\t" <<  temp_npoints << "\t" <<  temp_begintime << "\t" <<  temp_endtime << std::endl;
      if ( temp_feech<0 || temp_feech>255 )
      {
        std::cout << "ERROR, invalid FEECH " << temp_feech << " in MBD waveforms calibration" << std::endl;
        _status = -2;
        return _status;
      }

      _sherr_npts[temp_feech] = temp_npoints;
      _sherr_minrange[temp_feech] = temp_begintime;
      _sherr_maxrange[temp_feech] = temp_endtime;

      Double_t temp_val{0.};
      for (int isamp=0; isamp<temp_npoints; isamp++)
      {
        infile >> temp_val;
        _sherr_yerr[temp_feech].push_back( temp_val );
        if ( Verbosity() && (temp_feech==8 || temp_feech==255) )
        {
          std::cout << _sherr_yerr[temp_feech][isamp] << " ";
          if ( isamp%10==9 ) std::cout << std::endl;
        }
      }
      if ( Verbosity() ) std::cout << std::endl;
    }

    infile.close();
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;  // file not found
  }

  Verbosity(0);
  return 1;
}


int MbdCalib::Write_CDB_SampMax(const std::string& dbfile)
{
  CDBTTree *cdbttree {nullptr};

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version",1);
  cdbttree->CommitSingle();

  std::cout << "SAMPMAX" << std::endl;
  for (size_t ifeech=0; ifeech<_sampmax.size(); ifeech++) 
  {
    // store in a CDBTree
    cdbttree->SetIntValue(ifeech,"sampmax",_sampmax[ifeech]);

    if ( ifeech<5 || ifeech>=MbdDefs::MBD_N_FEECH-5 )
    {
      std::cout << ifeech << "\t" <<  cdbttree->GetIntValue(ifeech,"sampmax") << std::endl;
    }
  }

  cdbttree->Commit();
  //cdbttree->Print();

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}

int MbdCalib::Write_CDB_Shapes(const std::string& dbfile)
{
  CDBTTree *cdbttree {nullptr};

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version",1);
  cdbttree->CommitSingle();

  std::cout << "SHAPES" << std::endl;
  for (size_t ifeech=0; ifeech<_sampmax.size(); ifeech++) 
  {
    // store in a CDBTree
    cdbttree->SetIntValue(ifeech,"sampmax",_sampmax[ifeech]);

    if ( ifeech<5 || ifeech>=MbdDefs::MBD_N_FEECH-5 )
    {
      std::cout << ifeech << "\t" <<  cdbttree->GetIntValue(ifeech,"shape_npt") << std::endl;
    }
  }

  cdbttree->Commit();
  //cdbttree->Print();

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}

int MbdCalib::Write_CDB_All()
{
  return 1;
}

void MbdCalib::Reset()
{
  // Set all initial values
  _qfit_integ.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_mpv.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_sigma.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_integerr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_mpverr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_sigmaerr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_chi2ndf.fill(std::numeric_limits<float>::quiet_NaN());

  _tqfit_t0mean.fill(std::numeric_limits<float>::quiet_NaN());
  _tqfit_t0meanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _tqfit_t0sigma.fill(std::numeric_limits<float>::quiet_NaN());
  _tqfit_t0sigmaerr.fill(std::numeric_limits<float>::quiet_NaN());

  _sampmax.fill(-1);
}

