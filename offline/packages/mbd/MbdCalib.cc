#include "MbdCalib.h"

#include "MbdDefs.h"
#include "MbdGeomV1.h"

// Database Includes
#ifndef ONLINE
#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>
#endif

#include <phool/phool.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

#include <TString.h>

MbdCalib::MbdCalib()
{
  Reset();
  _mbdgeom = std::make_unique<MbdGeomV1>();

#ifndef ONLINE
  _rc = recoConsts::instance();
  if ( _rc->FlagExist("MBD_TEMPLATEFIT") )
  {
    do_templatefit = _rc->get_IntFlag("MBD_TEMPLATEFIT");
  }
  else
  {
    do_templatefit = 1;
  }
#else
  do_templatefit = 0;
#endif

}

int MbdCalib::Download_All()
{
  //std::cout << PHWHERE << " In MbdCalib::Download_All()" << std::endl;
  _status = 0;

  std::string bbc_caldir;

#ifdef ONLINE
  bbc_caldir = getenv("BBCCALIB");
  if (bbc_caldir.size()==0)
  {
    std::cout << "BBCCALIB environment variable not set" << std::endl;
    exit(1);
  }
#else
  if (Verbosity() > 0)
  {
    std::cout << "MBD CDB " << _rc->get_StringFlag("CDB_GLOBALTAG") << "\t" << _rc->get_uint64Flag("TIMESTAMP") << std::endl;
  }
  _cdb = CDBInterface::instance();

  if (_rc->FlagExist("MBD_CALDIR"))
  {
    bbc_caldir = _rc->get_StringFlag("MBD_CALDIR");
    std::cout << "Reading MBD Calibrations from " << bbc_caldir << std::endl;
  }

  // if rc flag MBD_CALDIR does not exist, we create it and set it to an empty string
  if (!_rc->FlagExist("MBD_CALDIR"))
  {
    Verbosity(0);
    std::string sampmax_url = _cdb->getUrl("MBD_SAMPMAX");
    if (Verbosity() > 0)
    {
      std::cout << "sampmax_url " << sampmax_url << std::endl;
    }
    Download_SampMax(sampmax_url);

    std::string qfit_url = _cdb->getUrl("MBD_QFIT");
    if (Verbosity() > 0)
    {
      std::cout << "qfit_url " << qfit_url << std::endl;
    }
    Download_Gains(qfit_url);

    std::string tt_t0_url = _cdb->getUrl("MBD_TT_T0");
    if ( Verbosity() > 0 )
    {
      std::cout << "tt_t0_url " << tt_t0_url << std::endl;
    }
    Download_TTT0(tt_t0_url);

    std::string tq_t0_url = _cdb->getUrl("MBD_TQ_T0");
    if (Verbosity() > 0)
    {
      std::cout << "tq_t0_url " << tq_t0_url << std::endl;
    }
    Download_TQT0(tq_t0_url);

    std::string timecorr_url = _cdb->getUrl("MBD_TIMECORR");
    if ( Verbosity() > 0 )
    {
      std::cout << "timecorr_url " << timecorr_url << std::endl;
    }
    Download_TimeCorr(timecorr_url);

    std::string slew_url = _cdb->getUrl("MBD_SLEWCORR");
    if ( Verbosity() > 0 )
    {
      std::cout << "slew_url " << slew_url << std::endl;
    }
    Download_SlewCorr(slew_url);

    if (do_templatefit)
    {
      std::string shape_url = _cdb->getUrl("MBD_SHAPES");
      if (Verbosity() > 0)
      {
        std::cout << "shape_url " << shape_url << std::endl;
      }
      Download_Shapes(shape_url);
    }
    Verbosity(0);
  }
#endif

  if ( bbc_caldir.size() > 0 )
  {
    std::string sampmax_file = bbc_caldir + "/mbd_sampmax.calib";
    Download_SampMax(sampmax_file);

    std::string qfit_file = bbc_caldir + "/mbd_qfit.calib";
    Download_Gains(qfit_file);

    std::string tq_t0_file = bbc_caldir + "/mbd_tq_t0.calib";
    Download_TQT0(tq_t0_file);

    std::string tt_tcorr_file = bbc_caldir + "/mbd_timecorr.calib";
    Download_TimeCorr(tt_tcorr_file);

    std::string tt_t0_file = bbc_caldir + "/mbd_tt_t0.calib";
    Download_TTT0(tt_t0_file);

    std::string slew_file = bbc_caldir + "/mbd_slewcorr.calib";
    Download_SlewCorr(slew_file);

    if (do_templatefit)
    {
      std::string shape_file = bbc_caldir + "/mbd_shape.calib";
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
  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from database
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
#endif

  if (dbase_file.EndsWith(".calib"))  // read from text file
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
  
  if ( std::isnan(_qfit_mpv[0]) )
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
  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from database
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
#endif

  if (dbase_file.EndsWith(".calib"))  // read from text file
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

  if ( std::isnan(_tqfit_t0mean[0]) )
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;
  }

  return 1;
}

int MbdCalib::Download_TTT0(const std::string& dbase_location)
{
  // Reset All Values
  _ttfit_t0mean.fill(std::numeric_limits<float>::quiet_NaN());
  _ttfit_t0meanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _ttfit_t0sigma.fill(std::numeric_limits<float>::quiet_NaN());
  _ttfit_t0sigmaerr.fill(std::numeric_limits<float>::quiet_NaN());

  if (Verbosity() > 0)
  {
    std::cout << "Opening " << dbase_location << std::endl;
  }
  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from database
  {
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    cdbttree->LoadCalibrations();

    for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
    {
      _ttfit_t0mean[ipmt] = cdbttree->GetFloatValue(ipmt, "ttfit_t0mean");
      _ttfit_t0meanerr[ipmt] = cdbttree->GetFloatValue(ipmt, "ttfit_t0meanerr");
      _ttfit_t0sigma[ipmt] = cdbttree->GetFloatValue(ipmt, "ttfit_t0sigma");
      _ttfit_t0sigmaerr[ipmt] = cdbttree->GetFloatValue(ipmt, "ttfit_t0sigmaerr");

      if (Verbosity() > 0)
      {
        if (ipmt < 5 || ipmt >= MbdDefs::MBD_N_PMT - 5)
        {
          std::cout << ipmt << "\t" << _ttfit_t0mean[ipmt] << std::endl;
        }
      }
    }
    delete cdbttree;
  }
#endif

  if (dbase_file.EndsWith(".calib"))  // read from text file
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
      infile >> _ttfit_t0mean[pmt] >> _ttfit_t0meanerr[pmt] >> _ttfit_t0sigma[pmt] >> _ttfit_t0sigmaerr[pmt];

      if (Verbosity() > 0)
      {
        if (pmt < 5 || pmt >= MbdDefs::MBD_N_PMT - 5)
        {
          std::cout << pmt << "\t" << _ttfit_t0mean[pmt] << "\t" << _ttfit_t0meanerr[pmt]
                    << "\t" << _ttfit_t0sigma[pmt] << "\t" << _ttfit_t0sigmaerr[pmt] << std::endl;
        }
      }
    }
    infile.close();
  }

  if ( std::isnan(_ttfit_t0mean[0]) )
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

  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from database
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
#endif

  if (dbase_file.EndsWith(".calib"))  // read from text file
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
          std::cout << "sampmax\t" << feech << "\t" << _sampmax[feech] << std::endl;
        }
      }
    }
    infile.close();
  }
  

  if ( _sampmax[0] == -1 )
  {
    std::cout << PHWHERE << ", ERROR, calib file not processed, " << dbase_location << std::endl;
    _status = -1;
    return _status;  // file not found
  }

  return 1;
}

int MbdCalib::Download_Shapes(const std::string& dbase_location)
{
  // Verbosity(100);
  if (Verbosity())
  {
    std::cout << "In MbdCalib::Download_Shapes" << std::endl;
  }
  // Reset All Values
  for (auto& shape : _shape_y)
  {
    shape.clear();
  }
  for (auto& sherr : _sherr_yerr)
  {
    sherr.clear();
  }

  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from database
  {
    if (Verbosity())
    {
      std::cout << "Reading from CDB " << dbase_location << std::endl;
    }
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    cdbttree->LoadCalibrations();

    for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
    {
      if (_mbdgeom->get_type(ifeech) == 0)
      {
        continue;  // skip t-channels
      }

      _shape_npts[ifeech] = cdbttree->GetIntValue(ifeech, "shape_npts");
      _shape_minrange[ifeech] = cdbttree->GetFloatValue(ifeech, "shape_min");
      _shape_maxrange[ifeech] = cdbttree->GetFloatValue(ifeech, "shape_max");

      _sherr_npts[ifeech] = cdbttree->GetIntValue(ifeech, "sherr_npts");
      _sherr_minrange[ifeech] = cdbttree->GetFloatValue(ifeech, "sherr_min");
      _sherr_maxrange[ifeech] = cdbttree->GetFloatValue(ifeech, "sherr_max");

      for (int ipt = 0; ipt < _shape_npts[ifeech]; ipt++)
      {
        int chtemp = 1000 * ipt + ifeech;

        float val = cdbttree->GetFloatValue(chtemp, "shape_val");
        _shape_y[ifeech].push_back(val);

        val = cdbttree->GetFloatValue(chtemp, "sherr_val");
        _sherr_yerr[ifeech].push_back(val);
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
#endif

  if (dbase_file.EndsWith(".calib"))  // read from text file
  {
    if (Verbosity())
    {
      std::cout << "Reading from " << dbase_location << std::endl;
    }
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
    while (infile >> temp_feech >> temp_npoints >> temp_begintime >> temp_endtime)
    {
      if (Verbosity())
      {
        std::cout << "shape " << temp_feech << "\t" << temp_npoints << "\t" << temp_begintime << "\t" << temp_endtime << std::endl;
      }
      if (temp_feech < 0 || temp_feech > 255)
      {
        std::cout << "ERROR, invalid FEECH " << temp_feech << " in MBD waveforms calibration" << std::endl;
        _status = -2;
        return _status;
      }

      _shape_npts[temp_feech] = temp_npoints;
      _shape_minrange[temp_feech] = temp_begintime;
      _shape_maxrange[temp_feech] = temp_endtime;

      float temp_val{0.};
      for (int isamp = 0; isamp < temp_npoints; isamp++)
      {
        infile >> temp_val;
        _shape_y[temp_feech].push_back(temp_val);
        if (Verbosity() && (temp_feech == 8 || temp_feech == 255))
        {
          std::cout << _shape_y[temp_feech][isamp] << " ";
          if (isamp % 10 == 9)
          {
            std::cout << std::endl;
          }
        }
      }
      if (Verbosity())
      {
        std::cout << std::endl;
      }
    }

    infile.close();

    // Now read in the sherr file
    std::string sherr_dbase_location = std::regex_replace(dbase_location, std::regex("bbc_shape.calib"), "bbc_sherr.calib");
    if (Verbosity())
    {
      std::cout << "Reading from " << sherr_dbase_location << std::endl;
    }
    infile.open(sherr_dbase_location);
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
    while (infile >> temp_feech >> temp_npoints >> temp_begintime >> temp_endtime)
    {
      if (Verbosity())
      {
        std::cout << "sheer " << temp_feech << "\t" << temp_npoints << "\t" << temp_begintime << "\t" << temp_endtime << std::endl;
      }
      if (temp_feech < 0 || temp_feech > 255)
      {
        std::cout << "ERROR, invalid FEECH " << temp_feech << " in MBD waveforms calibration" << std::endl;
        _status = -2;
        return _status;
      }

      _sherr_npts[temp_feech] = temp_npoints;
      _sherr_minrange[temp_feech] = temp_begintime;
      _sherr_maxrange[temp_feech] = temp_endtime;

      float temp_val{0.};
      for (int isamp = 0; isamp < temp_npoints; isamp++)
      {
        infile >> temp_val;
        _sherr_yerr[temp_feech].push_back(temp_val);
        if (Verbosity() && (temp_feech == 8 || temp_feech == 255))
        {
          std::cout << _sherr_yerr[temp_feech][isamp] << " ";
          if (isamp % 10 == 9)
          {
            std::cout << std::endl;
          }
        }
      }
      if (Verbosity())
      {
        std::cout << std::endl;
      }
    }

    infile.close();
  }

  if ( _shape_y[8].size()==0 )
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;  // file not found
  }

  // Verbosity(0);
  return 1;
}


int MbdCalib::Download_TimeCorr(const std::string& dbase_location)
{
  //Verbosity(100);
  if ( Verbosity() )
  {
    std::cout << "In MbdCalib::Download_TimeCorr" << std::endl;
  }
  // Reset All Values
  for(auto& tcorr : _tcorr_y) {
    tcorr.clear();
  }

  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from CDB database file
  {
    if ( Verbosity() )
    {
      std::cout << "Reading from CDB " << dbase_location << std::endl;
    }
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    cdbttree->LoadCalibrations();

    for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
    {
      if ( _mbdgeom->get_type(ifeech) == 1 )
      {
        continue;  // skip q-channels
      }

      _tcorr_npts[ifeech] = cdbttree->GetIntValue(ifeech, "tcorr_npts");
      _tcorr_minrange[ifeech] = cdbttree->GetFloatValue(ifeech, "tcorr_min");
      _tcorr_maxrange[ifeech] = cdbttree->GetFloatValue(ifeech, "tcorr_max");

      for (int ipt=0; ipt<_tcorr_npts[ifeech]; ipt++)
      {
        int chtemp = 1000*ipt + ifeech; // in cdbtree, entry has id = 1000*datapoint + ifeech

        float val = cdbttree->GetFloatValue(chtemp, "tcorr_val");
        _tcorr_y[ifeech].push_back( val );
      }

      if (Verbosity() > 0)
      {
        if (ifeech < 5 || ifeech >= MbdDefs::MBD_N_FEECH - 5)
        {
          std::cout << ifeech << "\t" << _tcorr_y[ifeech][0] << std::endl;
        }
      }
    }
    delete cdbttree;
  }
#endif

  if (dbase_file.EndsWith(".calib"))  // read from text file
  {
    if ( Verbosity() )
    {
      std::cout << "Reading from " << dbase_location << std::endl;
    }
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
      if ( Verbosity() )
      {
        std::cout << "tcorr " << temp_feech << "\t" <<  temp_npoints << "\t" <<  temp_begintime << "\t" <<  temp_endtime << std::endl;
      }
      if ( temp_feech<0 || temp_feech>255 )
      {
        std::cout << "ERROR, invalid FEECH " << temp_feech << " in MBD timecorr calibration" << std::endl;
        _status = -2;
        return _status;
      }

      _tcorr_npts[temp_feech] = temp_npoints;
      _tcorr_minrange[temp_feech] = temp_begintime;
      _tcorr_maxrange[temp_feech] = temp_endtime;

      float temp_val{0.};
      for (int isamp=0; isamp<temp_npoints; isamp++)
      {
        infile >> temp_val;
        _tcorr_y[temp_feech].push_back( temp_val );
        if ( Verbosity() && (temp_feech==0 || temp_feech==64) )
        {
          std::cout << _tcorr_y[temp_feech][isamp] << " ";
          if ( isamp%10==9 )
          {
            std::cout << std::endl;
          }
        }
      }
      if ( Verbosity() )
      {
        std::cout << std::endl;
      }
    }

    infile.close();
  }

  if ( _tcorr_y[0].size()==0 )
  {
    std::cout << PHWHERE << ", ERROR, MBD tcorr not loaded, " << dbase_location << std::endl;
    _status = -1;
    return _status;  // file not loaded
  }

  // Now we interpolate the timecorr
  for (size_t ifeech=0; ifeech<MbdDefs::MBD_N_FEECH; ifeech++) 
  {
    if ( _mbdgeom->get_type(ifeech) == 1 )
    {
      continue;  // skip q-channels
    }

    int step = static_cast<int>( (_tcorr_maxrange[ifeech] - _tcorr_minrange[ifeech]) / (_tcorr_npts[ifeech]-1) );
    //std::cout << ifeech << " step = " << step << std::endl;

    for (int itdc=0; itdc<=_tcorr_maxrange[ifeech]; itdc++)
    {
      int calib_index = itdc/step;
      int interp = itdc%step;

      // simple linear interpolation for now
      double slope = (_tcorr_y[ifeech][calib_index+1] - _tcorr_y[ifeech][calib_index])/step;
      float tcorr_interp = _tcorr_y[ifeech][calib_index] + interp*slope;
 
      _tcorr_y_interp[ifeech].push_back( tcorr_interp );

      if ( ifeech==0 && itdc<2*step && Verbosity() )
      {
        std::cout << _tcorr_y_interp[ifeech][itdc] << " ";
        if ( itdc%step==(step-1) )
        {
          std::cout << std::endl;
        }
      }
    }

  }

  //Verbosity(0);
  return 1;
}

int MbdCalib::Download_SlewCorr(const std::string& dbase_location)
{
  //Verbosity(100);
  if ( Verbosity() )
  {
    std::cout << "In MbdCalib::Download_SlewCorr" << std::endl;
  }
  // Reset All Values
  for(auto& scorr : _scorr_y) {
    scorr.clear();
  }
  std::fill(_scorr_npts.begin(), _scorr_npts.end(), 0);
  
  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from CDB database file
  {
    if ( Verbosity() )
    {
      std::cout << "Reading from CDB " << dbase_location << std::endl;
    }
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    cdbttree->LoadCalibrations();

    for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
    {
      if ( _mbdgeom->get_type(ifeech) == 1 )
      {
        continue;  // skip q-channels
      }

      _scorr_npts[ifeech] = cdbttree->GetIntValue(ifeech, "scorr_npts");
      _scorr_minrange[ifeech] = cdbttree->GetFloatValue(ifeech, "scorr_min");
      _scorr_maxrange[ifeech] = cdbttree->GetFloatValue(ifeech, "scorr_max");

      for (int ipt=0; ipt<_scorr_npts[ifeech]; ipt++)
      {
        int chtemp = 1000*ipt + ifeech; // in cdbtree, entry has id = 1000*datapoint + ifeech

        float val = cdbttree->GetFloatValue(chtemp, "scorr_val");
        _scorr_y[ifeech].push_back( val );
      }

      if (Verbosity() > 0)
      {
        if (ifeech < 5 || ifeech >= MbdDefs::MBD_N_FEECH - 5)
        {
          std::cout << ifeech << "\t" << _scorr_y[ifeech][0] << std::endl;
        }
      }
    }
    delete cdbttree;
  }
#endif

  if (dbase_file.EndsWith(".calib"))  // read from text file
  {
    if ( Verbosity() )
    {
      std::cout << "Reading from " << dbase_location << std::endl;
    }

    std::ifstream infile(dbase_location);
    if (!infile.is_open())
    {
      std::cout << PHWHERE << "unable to open " << dbase_location << std::endl;
      _status = -3;
      return _status;
    }

    int temp_feech = -1;
    int temp_npoints = 0;
    float temp_beginadc = -1;
    float temp_endadc = -1;
    while ( infile >> temp_feech >> temp_npoints >> temp_beginadc >> temp_endadc )
    {
      if ( Verbosity() )
      {
        std::cout << "scorr " << temp_feech << "\t" <<  temp_npoints << "\t" <<  temp_beginadc << "\t" <<  temp_endadc << std::endl;
      }

      if ( temp_feech<0 || temp_feech>255 )
      {
        std::cout << "ERROR, invalid FEECH " << temp_feech << " in MBD slewcorr calibration" << std::endl;
        _status = -2;
        return _status;
      }

      _scorr_npts[temp_feech] = temp_npoints;
      _scorr_minrange[temp_feech] = temp_beginadc;
      _scorr_maxrange[temp_feech] = temp_endadc;

      float temp_val{0.};
      for (int isamp=0; isamp<temp_npoints; isamp++)
      {
        infile >> temp_val;
        _scorr_y[temp_feech].push_back( temp_val );
        if ( Verbosity() && (temp_feech==0 || temp_feech==64) )
        {
          std::cout << _scorr_y[temp_feech][isamp] << " ";
          if ( isamp%10==9 )
          {
            std::cout << std::endl;
          }
        }
      }
      if ( Verbosity() )
      {
        std::cout << std::endl;
      }
    }

    infile.close();
  }

  if ( _scorr_y[0].size()==0 )
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;  // file not found
  }

  // Now we interpolate the slewcorr
  for (size_t ifeech=0; ifeech<MbdDefs::MBD_N_FEECH; ifeech++) 
  {
    if ( _mbdgeom->get_type(ifeech) == 1 )
    {
      continue;  // skip q-channels
    }
    // skip bad t-channels
    if ( _scorr_npts[ifeech] == 0 )
    {
      //std::cout << "skipping " << ifeech << std::endl;
      continue;
    }

    int step = static_cast<int>( (_scorr_maxrange[ifeech] - _scorr_minrange[ifeech]) / (_scorr_npts[ifeech]-1) );
    //std::cout << ifeech << " step = " << step << std::endl;

    for (int iadc=0; iadc<=_scorr_maxrange[ifeech]; iadc++)
    {
      int calib_index = iadc/step;
      int interp = iadc%step;

      // simple linear interpolation for now
      double slope = (_scorr_y[ifeech][calib_index+1] - _scorr_y[ifeech][calib_index])/step;
      float scorr_interp = _scorr_y[ifeech][calib_index] + interp*slope;
 
      _scorr_y_interp[ifeech].push_back( scorr_interp );


      if ( ifeech==4 && iadc<12 && Verbosity() )
      {
        if ( iadc==0 )
        {
          std::cout << "slewcorr " << ifeech << "\t" << _scorr_npts[ifeech] << "\t"
            << _scorr_minrange[ifeech] << "\t" << _scorr_maxrange[ifeech] << std::endl;
        }
        std::cout << _scorr_y_interp[ifeech][iadc] << " ";
        if ( iadc%step==(step-1) )
        {
          std::cout << std::endl;
        }
      }
    }

  }

  //Verbosity(0);
  return 1;
}

#ifndef ONLINE
int MbdCalib::Write_CDB_SampMax(const std::string& dbfile)
{
  CDBTTree* cdbttree{ nullptr };

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version", 1);
  cdbttree->CommitSingle();

  std::cout << "SAMPMAX" << std::endl;
  for (size_t ifeech = 0; ifeech < _sampmax.size(); ifeech++)
  {
    // store in a CDBTree
    cdbttree->SetIntValue(ifeech, "sampmax", _sampmax[ifeech]);

    if (ifeech < 12 || ifeech >= MbdDefs::MBD_N_FEECH - 5)
    {
      std::cout << ifeech << "\t" << cdbttree->GetIntValue(ifeech, "sampmax") << std::endl;
    }
  }

  cdbttree->Commit();
  // cdbttree->Print();

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

int MbdCalib::Write_SampMax(const std::string& dbfile)
{
  std::ofstream cal_file;
  cal_file.open(dbfile);
  for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
  {
    cal_file << ifeech << "\t" << _sampmax[ifeech] << std::endl;
  }
  cal_file.close();

  return 1;
}

#ifndef ONLINE
int MbdCalib::Write_CDB_TTT0(const std::string& dbfile)
{
  CDBTTree* cdbttree{ nullptr };

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version", 1);
  cdbttree->CommitSingle();

  std::cout << "TTT0" << std::endl;
  for (size_t ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    // store in a CDBTree
    cdbttree->SetFloatValue(ipmt, "ttfit_t0mean", _ttfit_t0mean[ipmt]);
    cdbttree->SetFloatValue(ipmt, "ttfit_t0meanerr", _ttfit_t0meanerr[ipmt]);
    cdbttree->SetFloatValue(ipmt, "ttfit_t0sigma", _ttfit_t0sigma[ipmt]);
    cdbttree->SetFloatValue(ipmt, "ttfit_t0sigmaerr", _ttfit_t0sigmaerr[ipmt]);

    if (ipmt < 5 || ipmt >= MbdDefs::MBD_N_PMT - 5)
    {
      std::cout << ipmt << "\t" << cdbttree->GetFloatValue(ipmt, "ttfit_t0mean") << std::endl;
    }
  }

  cdbttree->Commit();
  // cdbttree->Print();

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

int MbdCalib::Write_TTT0(const std::string& dbfile)
{
  std::ofstream cal_t0_file;
  cal_t0_file.open(dbfile);
  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    cal_t0_file << ipmt << "\t" << _ttfit_t0mean[ipmt] << "\t" << _ttfit_t0meanerr[ipmt]
      << "\t" << _ttfit_t0sigma[ipmt] << "\t" << _ttfit_t0sigmaerr[ipmt] << std::endl;
  }
  cal_t0_file.close();

  return 1;
}

#ifndef ONLINE
int MbdCalib::Write_CDB_TQT0(const std::string& dbfile)
{
  CDBTTree* cdbttree{ nullptr };

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version", 1);
  cdbttree->CommitSingle();

  std::cout << "TQT0" << std::endl;
  for (size_t ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    // store in a CDBTree
    cdbttree->SetFloatValue(ipmt, "tqfit_t0mean", _tqfit_t0mean[ipmt]);
    cdbttree->SetFloatValue(ipmt, "tqfit_t0meanerr", _tqfit_t0meanerr[ipmt]);
    cdbttree->SetFloatValue(ipmt, "tqfit_t0sigma", _tqfit_t0sigma[ipmt]);
    cdbttree->SetFloatValue(ipmt, "tqfit_t0sigmaerr", _tqfit_t0sigmaerr[ipmt]);

    if (ipmt < 5 || ipmt >= MbdDefs::MBD_N_PMT - 5)
    {
      std::cout << ipmt << "\t" << cdbttree->GetFloatValue(ipmt, "tqfit_t0mean") << std::endl;
    }
  }

  cdbttree->Commit();
  // cdbttree->Print();

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

int MbdCalib::Write_TQT0(const std::string& dbfile)
{
  std::ofstream cal_t0_file;
  cal_t0_file.open(dbfile);
  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    cal_t0_file << ipmt << "\t" << _tqfit_t0mean[ipmt] << "\t" << _tqfit_t0meanerr[ipmt]
      << "\t" << _tqfit_t0sigma[ipmt] << "\t" << _tqfit_t0sigmaerr[ipmt] << std::endl;
  }
  cal_t0_file.close();

  return 1;
}

#ifndef ONLINE
int MbdCalib::Write_CDB_Shapes(const std::string& dbfile)
{
  // store in a CDBTree
  CDBTTree* cdbttree{nullptr};

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree(dbfile);
  cdbttree->SetSingleIntValue("version", 1);
  cdbttree->CommitSingle();

  std::cout << "SHAPES" << std::endl;
  for (unsigned int ifeech = 0; ifeech < _sampmax.size(); ifeech++)
  {
    if (_mbdgeom->get_type(ifeech) == 0)
    {
      continue;  // skip t-channels
    }

    cdbttree->SetIntValue(ifeech, "shape_npts", _shape_npts[ifeech]);
    cdbttree->SetFloatValue(ifeech, "shape_min", _shape_minrange[ifeech]);
    cdbttree->SetFloatValue(ifeech, "shape_max", _shape_maxrange[ifeech]);

    for (int ipt = 0; ipt < _shape_npts[ifeech]; ipt++)
    {
      int temp_ch = ipt * 1000 + ifeech;
      cdbttree->SetFloatValue(temp_ch, "shape_val", _shape_y[ifeech][ipt]);
    }

    cdbttree->SetIntValue(ifeech, "sherr_npts", _sherr_npts[ifeech]);
    cdbttree->SetFloatValue(ifeech, "sherr_min", _sherr_minrange[ifeech]);
    cdbttree->SetFloatValue(ifeech, "sherr_max", _sherr_maxrange[ifeech]);

    for (int ipt = 0; ipt < _shape_npts[ifeech]; ipt++)
    {
      int temp_ch = ipt * 1000 + ifeech;
      cdbttree->SetFloatValue(temp_ch, "sherr_val", _sherr_yerr[ifeech][ipt]);
    }
  }

  cdbttree->Commit();
  // cdbttree->Print();

  for (unsigned int ifeech = 0; ifeech < _sampmax.size(); ifeech++)
  {
    if (_mbdgeom->get_type(ifeech) == 0)
    {
      continue;  // skip t-channels
    }

    if (ifeech < 5 || ifeech >= MbdDefs::MBD_N_FEECH - 5)
    {
      std::cout << ifeech << "\t" << cdbttree->GetIntValue(ifeech, "shape_npts") << std::endl;
      for (int ipt = 0; ipt < 10; ipt++)
      {
        int temp_ch = ipt * 1000 + (int)ifeech;
        std::cout << cdbttree->GetFloatValue(temp_ch, "shape_val") << "  ";
      }
      std::cout << std::endl;
    }
  }

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

#ifndef ONLINE
int MbdCalib::Write_CDB_TimeCorr(const std::string& dbfile)
{
  // store in a CDBTree
  CDBTTree *cdbttree {nullptr};

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version",1);
  cdbttree->CommitSingle();

  std::cout << "TIMECORR" << std::endl;
  for (int ifeech=0; ifeech<MbdDefs::MBD_N_FEECH; ifeech++) 
  {
    if ( _mbdgeom->get_type(ifeech) == 1 )
    {
      continue;  // skip q-channels
    }

    cdbttree->SetIntValue(ifeech,"tcorr_npts",_tcorr_npts[ifeech]);
    cdbttree->SetFloatValue(ifeech,"tcorr_min",_tcorr_minrange[ifeech]);
    cdbttree->SetFloatValue(ifeech,"tcorr_max",_tcorr_maxrange[ifeech]);

    for (int ipt=0; ipt<_tcorr_npts[ifeech]; ipt++)
    {
      int temp_ch = ipt*1000 + (int)ifeech;
      cdbttree->SetFloatValue(temp_ch,"tcorr_val",_tcorr_y[ifeech][ipt]);
    }
  }

  cdbttree->Commit();
  //cdbttree->Print();

  for (size_t ifeech=0; ifeech<MbdDefs::MBD_N_FEECH; ifeech++) 
  {
    if ( _mbdgeom->get_type(ifeech) == 1 )
    {
      continue;  // skip q-channels
    }

    if ( ifeech<5 || ifeech>=MbdDefs::MBD_N_FEECH-5-8 )
    {
      std::cout << ifeech << "\t" <<  cdbttree->GetIntValue(ifeech,"tcorr_npts") << std::endl;
      for (int ipt=0; ipt<10; ipt++)
      {
        int temp_ch = ipt*1000 + (int)ifeech;
        std::cout << cdbttree->GetFloatValue(temp_ch,"tcorr_val") << "  ";
      }
      std::cout << std::endl;
    }
  }

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

#ifndef ONLINE
int MbdCalib::Write_CDB_SlewCorr(const std::string& dbfile)
{
  // store in a CDBTree
  CDBTTree *cdbttree {nullptr};

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version",1);
  cdbttree->CommitSingle();

  std::cout << "SLEWCORR" << std::endl;
  //for (size_t ifeech=0; ifeech<_sampmax.size(); ifeech++) 
  for (size_t ifeech=0; ifeech<MbdDefs::MBD_N_FEECH; ifeech++) 
  {
    if ( _mbdgeom->get_type(ifeech) == 1 )
    {
      continue;  // skip q-channels
    }

    cdbttree->SetIntValue(ifeech,"scorr_npts",_scorr_npts[ifeech]);
    cdbttree->SetFloatValue(ifeech,"scorr_min",_scorr_minrange[ifeech]);
    cdbttree->SetFloatValue(ifeech,"scorr_max",_scorr_maxrange[ifeech]);

    for (int ipt=0; ipt<_scorr_npts[ifeech]; ipt++)
    {
      int temp_ch = ipt*1000 + (int)ifeech;
      cdbttree->SetFloatValue(temp_ch,"scorr_val",_scorr_y[ifeech][ipt]);
    }
  }

  cdbttree->Commit();
  //cdbttree->Print();

  for (size_t ifeech=0; ifeech<MbdDefs::MBD_N_FEECH; ifeech++) 
  {
    if ( _mbdgeom->get_type(ifeech) == 1 )
    {
      continue;  // skip q-channels
    }

    if ( ifeech<5 || ifeech>=MbdDefs::MBD_N_FEECH-8-5 )
    {
      std::cout << ifeech << "\t" <<  cdbttree->GetIntValue(ifeech,"scorr_npts") << std::endl;
      for (int ipt=0; ipt<10; ipt++)
      {
        int temp_ch = ipt*1000 + (int)ifeech;
        std::cout << cdbttree->GetFloatValue(temp_ch,"scorr_val") << "  ";
      }
      std::cout << std::endl;
    }
  }

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

#ifndef ONLINE
int MbdCalib::Write_CDB_Gains(const std::string& dbfile)
{
  CDBTTree* cdbttree{ nullptr };

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version", 1);
  cdbttree->CommitSingle();

  std::cout << "MBD_QFIT" << std::endl;
  for (size_t ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    // store in a CDBTree
    cdbttree->SetFloatValue(ipmt, "qfit_integ", _qfit_integ[ipmt]);
    cdbttree->SetFloatValue(ipmt, "qfit_mpv", _qfit_mpv[ipmt]);
    cdbttree->SetFloatValue(ipmt, "qfit_sigma", _qfit_sigma[ipmt]);
    cdbttree->SetFloatValue(ipmt, "qfit_integerr", _qfit_integerr[ipmt]);
    cdbttree->SetFloatValue(ipmt, "qfit_mpverr", _qfit_mpverr[ipmt]);
    cdbttree->SetFloatValue(ipmt, "qfit_sigmaerr", _qfit_sigmaerr[ipmt]);
    cdbttree->SetFloatValue(ipmt, "qfit_chi2ndf", _qfit_chi2ndf[ipmt]);

    if (ipmt < 5 || ipmt >= MbdDefs::MBD_N_PMT - 5)
    {
      std::cout << ipmt << "\t" << cdbttree->GetFloatValue(ipmt, "qfit_mpv") << std::endl;
    }
  }

  cdbttree->Commit();
  // cdbttree->Print();

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

#ifndef ONLINE
int MbdCalib::Write_CDB_All()
{
  return 1;
}
#endif

void MbdCalib::Update_TQT0(const float dz)
{
  // dz is what we need to move the MBD z by
  const float dt = dz/MbdDefs::C;

  for (int ipmt=0; ipmt<MbdDefs::MBD_N_PMT; ipmt++)
  {
    if ( ipmt<64 )  // south
    {
      _tqfit_t0mean[ipmt] -= dt;
    }
    else
    {
      _tqfit_t0mean[ipmt] += dt;
    }
  }
}

void MbdCalib::Update_TTT0(const float dz)
{
  // dz is what we need to move the MBD z by
  const float dt = dz/MbdDefs::C;

  for (int ipmt=0; ipmt<MbdDefs::MBD_N_PMT; ipmt++)
  {
    if ( ipmt<64 )  // south
    {
      _ttfit_t0mean[ipmt] -= dt;
    }
    else
    {
      _ttfit_t0mean[ipmt] += dt;
    }
  }
}

void MbdCalib::Reset_TTT0()
{
  _ttfit_t0mean.fill( 0. );
  _ttfit_t0meanerr.fill( 0. );
  _ttfit_t0sigma.fill(std::numeric_limits<float>::quiet_NaN());
  _ttfit_t0sigmaerr.fill(std::numeric_limits<float>::quiet_NaN());
}

void MbdCalib::Reset_TQT0()
{
  _tqfit_t0mean.fill( 0. );
  _tqfit_t0meanerr.fill( 0. );
  _tqfit_t0sigma.fill(std::numeric_limits<float>::quiet_NaN());
  _tqfit_t0sigmaerr.fill(std::numeric_limits<float>::quiet_NaN());
}

void MbdCalib::Reset_Gains()
{
  // Set all initial values
  _qfit_integ.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_mpv.fill( 1.0 );
  _qfit_sigma.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_integerr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_mpverr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_sigmaerr.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_chi2ndf.fill(std::numeric_limits<float>::quiet_NaN());
}

void MbdCalib::Reset()
{
  Reset_TTT0();
  Reset_TQT0();
  Reset_Gains();

  _sampmax.fill(-1);
}
