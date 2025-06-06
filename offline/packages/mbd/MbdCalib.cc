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
  // clang-tidy suggests to put this into the ctor but then ONLINE
  // does not compile anymore
  // NOLINTNEXTLINE(cppcoreguidelines-prefer-member-initializer)
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
  if ( Verbosity()>0 )
  {
    std::cout << PHWHERE << " In MbdCalib::Download_All()" << std::endl;
  }
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

    std::string t0corr_url = _cdb->getUrl("MBD_T0CORR");
    if ( Verbosity() > 0 )
    {
      std::cout << "t0corr_url " << t0corr_url << std::endl;
    }
    Download_T0Corr(t0corr_url);

    std::string ped_url = _cdb->getUrl("MBD_PED");
    if (Verbosity() > 0)
    {
      std::cout << "ped_url " << ped_url << std::endl;
    }
    Download_Ped(ped_url);

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

    std::string pileup_url = _cdb->getUrl("MBD_PILEUP");
    if (Verbosity() > 0)
    {
      std::cout << "pileup_url " << pileup_url << std::endl;
    }
    Download_Pileup(pileup_url);

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

  // download local calibs (text versions)
  if ( !bbc_caldir.empty() )
  {
    std::string sampmax_file = bbc_caldir + "/mbd_sampmax.calib";
    Download_SampMax(sampmax_file);

    std::string qfit_file = bbc_caldir + "/mbd_qfit.calib";
    Download_Gains(qfit_file);

    std::string tq_t0_file = bbc_caldir + "/mbd_tq_t0.calib";
    Download_TQT0(tq_t0_file);

    std::string tt_t0_file = bbc_caldir + "/mbd_tt_t0.calib";
    Download_TTT0(tt_t0_file);

    std::string t0corr_file = bbc_caldir + "/mbd_t0corr.calib";
    Download_T0Corr(t0corr_file);

    std::string ped_file = bbc_caldir + "/mbd_ped.calib";
    Download_Ped(ped_file);

    std::string tt_tcorr_file = bbc_caldir + "/mbd_timecorr.calib";
    Download_TimeCorr(tt_tcorr_file);

    std::string slew_file = bbc_caldir + "/mbd_slewcorr.calib";
    Download_SlewCorr(slew_file);

    std::string pileup_file = bbc_caldir + "/mbd_pileup.calib";
    Download_Pileup(pileup_file);

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


int MbdCalib::Download_T0Corr(const std::string& dbase_location)
{
  // Reset All Values
  _t0corrmean = 0.;
  _t0corrmeanerr = 0.;
  _t0corr_fitmean.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_fitmeanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_fitsigma.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_fitsigmaerr.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_hmean.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_hmeanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_hstddev.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_hstddeverr.fill(std::numeric_limits<float>::quiet_NaN());

  if ( dbase_location.empty() )
  {
    return 0;
  }

  if (Verbosity() > 0)
  {
    std::cout << "Opening " << dbase_location << std::endl;
  }
  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from database
  {
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    if ( cdbttree == nullptr )
    {
      std::cerr << "T0Corr not found, skipping" << std::endl;
      _status = -1;
      return _status;
    }
    cdbttree->LoadCalibrations();

    _t0corrmean = cdbttree->GetSingleFloatValue("_t0corrmean");
    _t0corrmeanerr = cdbttree->GetSingleFloatValue("_t0corrmeanerr");
    _t0corr_fitmean[0] = cdbttree->GetSingleFloatValue("_t0corr_fitmean0");
    _t0corr_fitmeanerr[0] = cdbttree->GetSingleFloatValue("_t0corr_fitmeanerr0");
    _t0corr_fitsigma[0] = cdbttree->GetSingleFloatValue("_t0corr_fitsigma0");
    _t0corr_fitsigmaerr[0] = cdbttree->GetSingleFloatValue("_t0corr_fitsigmaerr0");
    _t0corr_fitmean[1] = cdbttree->GetSingleFloatValue("_t0corr_fitmean1");
    _t0corr_fitmeanerr[1] = cdbttree->GetSingleFloatValue("_t0corr_fitmeanerr1");
    _t0corr_fitsigma[1] = cdbttree->GetSingleFloatValue("_t0corr_fitsigma1");
    _t0corr_fitsigmaerr[1] = cdbttree->GetSingleFloatValue("_t0corr_fitsigmaerr1");
    _t0corr_hmean[0] = cdbttree->GetSingleFloatValue("_t0corr_hmean0");
    _t0corr_hmeanerr[0] = cdbttree->GetSingleFloatValue("_t0corr_hmeanerr0");
    _t0corr_hstddev[0] = cdbttree->GetSingleFloatValue("_t0corr_hstddev0");
    _t0corr_hstddeverr[0] = cdbttree->GetSingleFloatValue("_t0corr_hstddeverr0");
    _t0corr_hmean[1] = cdbttree->GetSingleFloatValue("_t0corr_hmean1");
    _t0corr_hmeanerr[1] = cdbttree->GetSingleFloatValue("_t0corr_hmeanerr1");
    _t0corr_hstddev[1] = cdbttree->GetSingleFloatValue("_t0corr_hstddev1");
    _t0corr_hstddeverr[1] = cdbttree->GetSingleFloatValue("_t0corr_hstddeverr1");

    if (Verbosity() > 0)
    {
      std::cout << "T0Corr\t" << _t0corrmean << std::endl;
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

    infile >> _t0corrmean >> _t0corrmeanerr
      >> _t0corr_fitmean[0] >> _t0corr_fitmeanerr[0] >> _t0corr_fitsigma[0] >> _t0corr_fitsigmaerr[0]
      >> _t0corr_fitmean[1] >> _t0corr_fitmeanerr[1] >> _t0corr_fitsigma[1] >> _t0corr_fitsigmaerr[1]
      >> _t0corr_hmean[0] >> _t0corr_hmeanerr[0] >> _t0corr_hstddev[0] >> _t0corr_hstddeverr[0]
      >> _t0corr_hmean[1] >> _t0corr_hmeanerr[1] >> _t0corr_hstddev[1] >> _t0corr_hstddeverr[1];

    if (Verbosity() > 0)
    {
      std::cout << "T0Corr\t" << _t0corrmean << "\t" << _t0corrmeanerr << std::endl;
    }
    infile.close();
  }

  if ( std::isnan(_t0corrmean) )
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;
  }

  return 1;
}

int MbdCalib::Download_Ped(const std::string& dbase_location)
{
  // Reset All Values
  _pedmean.fill(std::numeric_limits<float>::quiet_NaN());
  _pedmeanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _pedsigma.fill(std::numeric_limits<float>::quiet_NaN());
  _pedsigmaerr.fill(std::numeric_limits<float>::quiet_NaN());

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

    for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
    {
      _pedmean[ifeech] = cdbttree->GetFloatValue(ifeech, "pedmean");
      _pedmeanerr[ifeech] = cdbttree->GetFloatValue(ifeech, "pedmeanerr");
      _pedsigma[ifeech] = cdbttree->GetFloatValue(ifeech, "pedsigma");
      _pedsigmaerr[ifeech] = cdbttree->GetFloatValue(ifeech, "pedsigmaerr");

      if (Verbosity() > 0)
      {
        if (ifeech < 5 || ifeech >= MbdDefs::MBD_N_FEECH - 5)
        {
          std::cout << ifeech << "\t" << _pedmean[ifeech] << std::endl;
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
      infile >> _pedmean[feech] >> _pedmeanerr[feech] >> _pedsigma[feech] >> _pedsigmaerr[feech];

      if (Verbosity() > 0)
      {
        if (feech < 5 || feech >= MbdDefs::MBD_N_FEECH - 5)
        {
          std::cout << feech << "\t" << _pedmean[feech] << "\t" << _pedmeanerr[feech]
                    << "\t" << _pedsigma[feech] << "\t" << _pedsigmaerr[feech] << std::endl;
        }
      }
    }
    infile.close();
  }

  if ( std::isnan(_pedmean[0]) )
  {
    std::cout << PHWHERE << ", WARNING, ped calib missing, " << dbase_location << std::endl;
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
    std::cout << PHWHERE << ", WARNING, sampmax calib missing, " << dbase_location << std::endl;
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
        int chtemp = (1000 * ipt) + ifeech;

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

  if ( _shape_y[8].empty() )
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
        int chtemp = (1000*ipt) + ifeech; // in cdbtree, entry has id = 1000*datapoint + ifeech

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

  if ( _tcorr_y[0].empty() )
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
      float tcorr_interp = _tcorr_y[ifeech][calib_index] + (interp*slope);
 
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
        int chtemp = (1000*ipt) + ifeech; // in cdbtree, entry has id = 1000*datapoint + ifeech

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

  if ( _scorr_y[0].empty() )
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
      float scorr_interp = _scorr_y[ifeech][calib_index] + (interp*slope);
 
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

int MbdCalib::Download_Pileup(const std::string& dbase_location)
{
  // Reset All Values
  Reset_Pileup();

  if (Verbosity() > 0)
  {
    std::cout << "Opening " << dbase_location << std::endl;
  }
  TString dbase_file = dbase_location;

#ifndef ONLINE
  if (dbase_file.EndsWith(".root"))  // read from database
  {
    CDBTTree* cdbttree = new CDBTTree(dbase_location);
    if ( cdbttree == nullptr )
    {
      std::cerr << "MBD pileup calib not found, skipping" << std::endl;
      _status = -1;
      return _status;
    }
    cdbttree->LoadCalibrations();

    for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
    {
      _pileup_p0[ifeech] = cdbttree->GetFloatValue(ifeech, "pileup_p0");
      _pileup_p0err[ifeech] = cdbttree->GetFloatValue(ifeech, "pileup_p0err");
      _pileup_p1[ifeech] = cdbttree->GetFloatValue(ifeech, "pileup_p1");
      _pileup_p1err[ifeech] = cdbttree->GetFloatValue(ifeech, "pileup_p1err");
      _pileup_p2[ifeech] = cdbttree->GetFloatValue(ifeech, "pileup_p2");
      _pileup_p2err[ifeech] = cdbttree->GetFloatValue(ifeech, "pileup_p2err");
      _pileup_chi2ndf[ifeech] = cdbttree->GetFloatValue(ifeech, "pileup_chi2ndf");
      if (Verbosity() > 0)
      {
        if (ifeech < 2 || ifeech >= (MbdDefs::MBD_N_FEECH-2) )
        {
          std::cout << ifeech << "\t" << _pileup_p0[ifeech] << std::endl;
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
      infile >> _pileup_p0[feech] >> _pileup_p0err[feech]
             >> _pileup_p1[feech] >> _pileup_p1err[feech]
             >> _pileup_p2[feech] >> _pileup_p2err[feech]
             >> _pileup_chi2ndf[feech];

      if (Verbosity() > 0)
      {
        if (feech < 2 || feech >= MbdDefs::MBD_N_PMT - 2)
        {
          std::cout << feech << "\t" << _pileup_p0[feech] << "\t" << _pileup_p0err[feech]
                             << "\t" << _pileup_p1[feech] << "\t" << _pileup_p1err[feech]
                             << "\t" << _pileup_p2[feech] << "\t" << _pileup_p2err[feech]
                             << "\t" << _pileup_chi2ndf[feech] << std::endl;
        }
      }
    }
  }
  
  if ( std::isnan(_pileup_p0[0]) )
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;
  }

  return 1;
}

int MbdCalib::Download_Thresholds(const std::string& dbase_location)
{
  // Reset All Values
  Reset_Thresholds();

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
      _thresh_mean[ipmt] = cdbttree->GetFloatValue(ipmt, "thresh_mean");
      _thresh_meanerr[ipmt] = cdbttree->GetFloatValue(ipmt, "thresh_meanerr");
      _thresh_width[ipmt] = cdbttree->GetFloatValue(ipmt, "thresh_width");
      _thresh_widtherr[ipmt] = cdbttree->GetFloatValue(ipmt, "thresh_widtherr");
      _thresh_eff[ipmt] = cdbttree->GetFloatValue(ipmt, "thresh_eff");
      _thresh_efferr[ipmt] = cdbttree->GetFloatValue(ipmt, "thresh_efferr");
      _thresh_chi2ndf[ipmt] = cdbttree->GetFloatValue(ipmt, "thresh_chi2ndf");
      if (Verbosity() > 0)
      {
        if (ipmt < 5)
        {
          std::cout << ipmt << "\t" << _thresh_mean[ipmt] << std::endl;
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
      infile >> _thresh_mean[pmt] >> _thresh_meanerr[pmt] >> _thresh_width[pmt] >> _thresh_widtherr[pmt] >> _thresh_eff[pmt] >> _thresh_efferr[pmt] >> _thresh_chi2ndf[pmt];
      if (Verbosity() > 0)
      {
        if (pmt < 5 || pmt >= MbdDefs::MBD_N_PMT - 5)
        {
          std::cout << pmt << "\t" << _thresh_mean[pmt] << "\t" << _thresh_meanerr[pmt] << "\t" << _thresh_width[pmt]
                    << "\t" << _thresh_widtherr[pmt] << "\t" << _thresh_eff[pmt] << "\t" << _thresh_efferr[pmt]
                    << "\t" << _thresh_chi2ndf[pmt] << std::endl;
        }
      }
    }
  }
  
  if ( std::isnan(_thresh_mean[0]) )
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << std::endl;
    _status = -1;
    return _status;
  }

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
int MbdCalib::Write_CDB_T0Corr(const std::string& dbfile)
{
  CDBTTree* cdbttree{ nullptr };

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version", 1);

  std::cout << "T0Corr" << std::endl;

  // store in a CDBTree
  cdbttree->SetSingleFloatValue("_t0corrmean", _t0corrmean);
  cdbttree->SetSingleFloatValue("_t0corrmeanerr", _t0corrmeanerr);
  cdbttree->SetSingleFloatValue("_t0corr_fitmean0", _t0corr_fitmean[0]);
  cdbttree->SetSingleFloatValue("_t0corr_fitmeanerr0", _t0corr_fitmeanerr[0]);
  cdbttree->SetSingleFloatValue("_t0corr_fitsigma0", _t0corr_fitsigma[0]);
  cdbttree->SetSingleFloatValue("_t0corr_fitsigmaerr0", _t0corr_fitsigmaerr[0]);
  cdbttree->SetSingleFloatValue("_t0corr_fitmean1", _t0corr_fitmean[1]);
  cdbttree->SetSingleFloatValue("_t0corr_fitmeanerr1", _t0corr_fitmeanerr[1]);
  cdbttree->SetSingleFloatValue("_t0corr_fitsigma1", _t0corr_fitsigma[1]);
  cdbttree->SetSingleFloatValue("_t0corr_fitsigmaerr1", _t0corr_fitsigmaerr[1]);
  cdbttree->SetSingleFloatValue("_t0corr_hmean0", _t0corr_hmean[0]);
  cdbttree->SetSingleFloatValue("_t0corr_hmeanerr0", _t0corr_hmeanerr[0]);
  cdbttree->SetSingleFloatValue("_t0corr_hstddev0", _t0corr_hstddev[0]);
  cdbttree->SetSingleFloatValue("_t0corr_hstddeverr0", _t0corr_hstddeverr[0]);
  cdbttree->SetSingleFloatValue("_t0corr_hmean1", _t0corr_hmean[1]);
  cdbttree->SetSingleFloatValue("_t0corr_hmeanerr1", _t0corr_hmeanerr[1]);
  cdbttree->SetSingleFloatValue("_t0corr_hstddev1", _t0corr_hstddev[1]);
  cdbttree->SetSingleFloatValue("_t0corr_hstddeverr1", _t0corr_hstddeverr[1]);

  std::cout << "T0Corr\t" << cdbttree->GetSingleFloatValue("_t0corrmean") << std::endl;

  cdbttree->CommitSingle();
  //cdbttree->Commit();
  //cdbttree->Print();

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

int MbdCalib::Write_T0Corr(const std::string& dbfile)
{
  std::ofstream cal_t0corr_file;
  cal_t0corr_file.open(dbfile);
  cal_t0corr_file << _t0corrmean << "\t" << _t0corrmeanerr << std::endl;
  cal_t0corr_file << _t0corr_fitmean[0] << "\t" << _t0corr_fitmeanerr[0] << "\t" << _t0corr_fitsigma[0] << "\t" << _t0corr_fitsigmaerr[0] << std::endl;
  cal_t0corr_file << _t0corr_fitmean[1] << "\t" << _t0corr_fitmeanerr[1] << "\t" << _t0corr_fitsigma[1] << "\t" << _t0corr_fitsigmaerr[1] << std::endl;
  cal_t0corr_file << _t0corr_hmean[0] << "\t" << _t0corr_hmeanerr[0] << "\t" << _t0corr_hstddev[0] << "\t" << _t0corr_hstddeverr[0] << std::endl;
  cal_t0corr_file << _t0corr_hmean[1] << "\t" << _t0corr_hmeanerr[1] << "\t" << _t0corr_hstddev[1] << "\t" << _t0corr_hstddeverr[1] << std::endl;
  cal_t0corr_file.close();

  return 1;
}

#ifndef ONLINE
int MbdCalib::Write_CDB_Ped(const std::string& dbfile)
{
  CDBTTree* cdbttree{ nullptr };

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version", 1);

  std::cout << "Ped" << std::endl;
  for (size_t ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
  {
    // store in a CDBTree
    cdbttree->SetFloatValue(ifeech, "pedmean", _pedmean[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pedmeanerr", _pedmeanerr[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pedsigma", _pedsigma[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pedsigmaerr", _pedsigmaerr[ifeech]);

    if (ifeech < 5 || ifeech >= MbdDefs::MBD_N_FEECH - 5)
    {
      std::cout << ifeech << "\t" << cdbttree->GetFloatValue(ifeech, "pedmean") << std::endl;
    }
  }
  cdbttree->CommitSingle();

  cdbttree->Commit();
  // cdbttree->Print();

  // for now we create the tree after reading it
  cdbttree->WriteCDBTTree();
  delete cdbttree;

  return 1;
}
#endif

int MbdCalib::Write_Ped(const std::string& dbfile)
{
  std::ofstream cal_ped_file;
  cal_ped_file.open(dbfile);
  for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
  {
    cal_ped_file << ifeech << "\t" << _pedmean[ifeech] << "\t" << _pedmeanerr[ifeech]
      << "\t" << _pedsigma[ifeech] << "\t" << _pedsigmaerr[ifeech] << std::endl;
  }
  cal_ped_file.close();

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
      int temp_ch = (ipt * 1000) + ifeech;
      cdbttree->SetFloatValue(temp_ch, "shape_val", _shape_y[ifeech][ipt]);
    }

    cdbttree->SetIntValue(ifeech, "sherr_npts", _sherr_npts[ifeech]);
    cdbttree->SetFloatValue(ifeech, "sherr_min", _sherr_minrange[ifeech]);
    cdbttree->SetFloatValue(ifeech, "sherr_max", _sherr_maxrange[ifeech]);

    for (int ipt = 0; ipt < _shape_npts[ifeech]; ipt++)
    {
      int temp_ch = (ipt * 1000) + ifeech;
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
        int temp_ch = (ipt * 1000) + (int)ifeech;
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
      int temp_ch = (ipt*1000) + ifeech;
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
        int temp_ch = (ipt*1000) + (int)ifeech;
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
      int temp_ch = (ipt*1000) + (int)ifeech;
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
        int temp_ch = (ipt*1000) + (int)ifeech;
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

int MbdCalib::Write_Gains(const std::string& dbfile)
{
  std::ofstream cal_gains_file;
  cal_gains_file.open(dbfile);
  for (int ipmtch = 0; ipmtch < MbdDefs::MBD_N_PMT; ipmtch++)
  {
    cal_gains_file << ipmtch << "\t" << _qfit_integ[ipmtch] << "\t" << _qfit_mpv[ipmtch]
      << "\t" << _qfit_sigma[ipmtch] << "\t" << _qfit_integerr[ipmtch]
      << "\t" << _qfit_mpverr[ipmtch] << "\t" << _qfit_sigmaerr[ipmtch]
      << "\t" << _qfit_chi2ndf[ipmtch] << std::endl;
  }
  cal_gains_file.close();

  return 1;
}

#ifndef ONLINE
int MbdCalib::Write_CDB_Pileup(const std::string& dbfile)
{
  CDBTTree* cdbttree{ nullptr };

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version", 1);
  cdbttree->CommitSingle();

  std::cout << "MBD_PILEUP" << std::endl;
  for (size_t ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
  {
    // store in a CDBTree
    cdbttree->SetFloatValue(ifeech, "pileup_p0", _pileup_p0[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pileup_p0err", _pileup_p0err[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pileup_p1", _pileup_p1[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pileup_p1err", _pileup_p1err[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pileup_p2", _pileup_p2[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pileup_p2err", _pileup_p2err[ifeech]);
    cdbttree->SetFloatValue(ifeech, "pileup_chi2ndf", _pileup_chi2ndf[ifeech]);

    if (ifeech < 5 || ifeech >= MbdDefs::MBD_N_PMT - 5)
    {
      std::cout << ifeech << "\t" << cdbttree->GetFloatValue(ifeech, "pileup_p0") << std::endl;
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

int MbdCalib::Write_Pileup(const std::string& dbfile)
{
  std::ofstream cal_pileup_file;
  cal_pileup_file.open(dbfile);
  for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
  {
    cal_pileup_file << ifeech << "\t" << _pileup_p0[ifeech] << "\t" << _pileup_p0err[ifeech]
      << "\t" << _pileup_p1[ifeech] << "\t" << _pileup_p1err[ifeech]
      << "\t" << _pileup_p2[ifeech] << "\t" << _pileup_p2err[ifeech]
      << "\t" << _pileup_chi2ndf[ifeech] << std::endl;
  }
  cal_pileup_file.close();

  return 1;
}

#ifndef ONLINE
int MbdCalib::Write_CDB_Thresholds(const std::string& dbfile)
{
  CDBTTree* cdbttree{ nullptr };

  std::cout << "Creating " << dbfile << std::endl;
  cdbttree = new CDBTTree( dbfile );
  cdbttree->SetSingleIntValue("version", 1);
  cdbttree->CommitSingle();

  std::cout << "MBD_THRESHOLDS" << std::endl;
  for (size_t ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    // store in a CDBTree
    cdbttree->SetFloatValue(ipmt, "thresh_mean", _thresh_mean[ipmt]);
    cdbttree->SetFloatValue(ipmt, "thresh_meanerr", _thresh_meanerr[ipmt]);
    cdbttree->SetFloatValue(ipmt, "thresh_width", _thresh_width[ipmt]);
    cdbttree->SetFloatValue(ipmt, "thresh_widtherr", _thresh_widtherr[ipmt]);
    cdbttree->SetFloatValue(ipmt, "thresh_eff", _thresh_eff[ipmt]);
    cdbttree->SetFloatValue(ipmt, "thresh_efferr", _thresh_efferr[ipmt]);
    cdbttree->SetFloatValue(ipmt, "thresh_chi2ndf", _thresh_chi2ndf[ipmt]);

    if (ipmt < 5 || ipmt >= MbdDefs::MBD_N_PMT - 5)
    {
      std::cout << ipmt << "\t" << cdbttree->GetFloatValue(ipmt, "thresh_mpv") << std::endl;
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

int MbdCalib::Write_Thresholds(const std::string& dbfile)
{
  std::ofstream cal_thresh_file;
  cal_thresh_file.open(dbfile);
  for (int ipmtch = 0; ipmtch < MbdDefs::MBD_N_PMT; ipmtch++)
  {
    cal_thresh_file << ipmtch << "\t" << _thresh_mean[ipmtch] << "\t" << _thresh_meanerr[ipmtch]
      << "\t" << _thresh_width[ipmtch] << "\t" << _thresh_widtherr[ipmtch]
      << "\t" << _thresh_eff[ipmtch] << "\t" << _thresh_efferr[ipmtch]
      << "\t" << _thresh_chi2ndf[ipmtch] << std::endl;
  }
  cal_thresh_file.close();

  return 1;
}

#ifndef ONLINE
int MbdCalib::Write_CDB_All()
{
  return 1;
}
#endif

// dz is what we need to move the MBD z by
// dt is what we change the MBD t0 by
void MbdCalib::Update_TQT0(const float dz, const float dt)
{
  const float z_dt = dz/MbdDefs::C;

  for (int ipmt=0; ipmt<MbdDefs::MBD_N_PMT; ipmt++)
  {
    // update zvtx
    if ( ipmt<64 )  // south
    {
      _tqfit_t0mean[ipmt] -= z_dt;
    }
    else
    {
      _tqfit_t0mean[ipmt] += z_dt;
    }

    // update t0
    _tqfit_t0mean[ipmt] -= dt;
  }
}

void MbdCalib::Update_TTT0(const float dz, const float dt)
{
  // dz is what we need to move the MBD z by
  const float z_dt = dz/MbdDefs::C;

  for (int ipmt=0; ipmt<MbdDefs::MBD_N_PMT; ipmt++)
  {
    if ( ipmt<64 )  // south
    {
      _ttfit_t0mean[ipmt] -= z_dt;
    }
    else
    {
      _ttfit_t0mean[ipmt] += z_dt;
    }

    // update t0
    _ttfit_t0mean[ipmt] -= dt;
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

void MbdCalib::Reset_T0Corr()
{
  _t0corrmean = 0.;
  _t0corrmeanerr = 0.;
  _t0corr_fitmean.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_fitmeanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_fitsigma.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_fitsigmaerr.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_hmean.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_hmeanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_hstddev.fill(std::numeric_limits<float>::quiet_NaN());
  _t0corr_hstddeverr.fill(std::numeric_limits<float>::quiet_NaN());
}

void MbdCalib::Reset_Ped()
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

void MbdCalib::Reset_Pileup()
{
  // Set all initial values
  _pileup_p0.fill(std::numeric_limits<float>::quiet_NaN());
  _pileup_p1.fill(std::numeric_limits<float>::quiet_NaN());
  _pileup_p2.fill(std::numeric_limits<float>::quiet_NaN());
  _pileup_p0err.fill(std::numeric_limits<float>::quiet_NaN());
  _pileup_p1err.fill(std::numeric_limits<float>::quiet_NaN());
  _pileup_p2err.fill(std::numeric_limits<float>::quiet_NaN());
  _qfit_chi2ndf.fill(std::numeric_limits<float>::quiet_NaN());
}

void MbdCalib::Reset_Thresholds()
{
  // Set all initial values
  _thresh_mean.fill(std::numeric_limits<float>::quiet_NaN());
  _thresh_meanerr.fill(std::numeric_limits<float>::quiet_NaN());
  _thresh_width.fill(std::numeric_limits<float>::quiet_NaN());
  _thresh_widtherr.fill(std::numeric_limits<float>::quiet_NaN());
  _thresh_eff.fill(std::numeric_limits<float>::quiet_NaN());
  _thresh_efferr.fill(std::numeric_limits<float>::quiet_NaN());
  _thresh_chi2ndf.fill(std::numeric_limits<float>::quiet_NaN());
}

void MbdCalib::Reset()
{
  Reset_TTT0();
  Reset_TQT0();
  Reset_Ped();
  Reset_Gains();
  Reset_T0Corr();
  Reset_Pileup();
  Reset_Thresholds();

  _sampmax.fill(-1);
}

void MbdCalib::set_ped(const int ifeech, const float m, const float merr, const float s, const float serr)
{
  _pedmean[ifeech] = m;
  _pedmeanerr[ifeech] = merr;
  _pedsigma[ifeech] = s;
  _pedsigmaerr[ifeech] = serr;
}


float MbdCalib::get_threshold(const int pmtch, const int rel_or_abs)
{
  if ( rel_or_abs==0 )
  {
    return _thresh_mean[pmtch]/_qfit_mpv[pmtch];
  }
  else
  {
    return _thresh_mean[pmtch];
  }

  return -1.; // error
}
