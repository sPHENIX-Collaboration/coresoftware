#include "BbcCalib.h"

// Database Includes
#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <phool/phool.h>
#include <cstring>

#include <TString.h>

using namespace std;
using namespace BbcDefs;

BbcCalib::BbcCalib()
{
  Reset();
  _rc = recoConsts::instance();
  _cdb = CDBInterface::instance();
}

int BbcCalib::Download_All()
{
  _status = 0;

  cout << "BBC CDB " << _rc->get_StringFlag("CDB_GLOBALTAG") << "\t" << _rc->get_uint64Flag("TIMESTAMP") << endl;

  _cdb = CDBInterface::instance();

  string bbc_caldir = _rc->get_StringFlag("BBC_CALDIR");
  if ( bbc_caldir.empty() )
  {
    string sampmax_url = _cdb->getUrl("BBC_SAMPMAX");
    cout << "sampmax_url " << sampmax_url << endl;
    Download_SampMax( sampmax_url );

    string qfit_url = _cdb->getUrl("BBC_QFIT");
    cout << "qfit_url " << qfit_url << endl;
    Download_Gains( qfit_url );

    string tq_t0_url = _cdb->getUrl("BBC_TQ_T0");
    cout << "tq_t0_url " << tq_t0_url << endl;
    Download_TQT0( tq_t0_url );
  }
  else
  {
    std::string sampmax_file = bbc_caldir + "/bbc_sampmax.calib";
    Download_SampMax( sampmax_file );
    std::string qfit_file = bbc_caldir + "/bbc_qfit.calib";
    Download_Gains( qfit_file );
    std::string tq_t0_file = bbc_caldir + "/bbc_tq_t0.calib";
    Download_TQT0( tq_t0_file );
  }

  return _status;
}

int BbcCalib::Download_Gains(const string& dbase_location)
{
  // Reset All Values
  for (int ipmt=0; ipmt<BBC_N_PMT; ipmt++)
  {
    _qfit_integ[ipmt] = NAN;
    _qfit_mpv[ipmt] = NAN;
    _qfit_sigma[ipmt] = NAN;
    _qfit_integerr[ipmt] = NAN;
    _qfit_mpverr[ipmt] = NAN;
    _qfit_sigmaerr[ipmt] = NAN;
    _qfit_chi2ndf[ipmt] = NAN;
  }

  cout << "Opening " << dbase_location << endl;
  TString dbase_file = dbase_location;
  if ( dbase_file.EndsWith(".root") )       // read from database
  {
    CDBTTree *cdbttree = new CDBTTree( dbase_location );
    cdbttree->LoadCalibrations();

    for (int ipmt=0; ipmt<BBC_N_PMT; ipmt++)
    {
      _qfit_integ[ipmt] = cdbttree->GetFloatValue(ipmt,"qfit_integ");
      _qfit_mpv[ipmt] = cdbttree->GetFloatValue(ipmt,"qfit_mpv");
      _qfit_sigma[ipmt] = cdbttree->GetFloatValue(ipmt,"qfit_sigma");
      _qfit_integerr[ipmt] = cdbttree->GetFloatValue(ipmt,"qfit_integerr");
      _qfit_mpverr[ipmt] = cdbttree->GetFloatValue(ipmt,"qfit_mpverr");
      _qfit_sigmaerr[ipmt] = cdbttree->GetFloatValue(ipmt,"qfit_sigmaerr");
      _qfit_chi2ndf[ipmt] = cdbttree->GetFloatValue(ipmt,"qfit_chi2ndf");
      if (ipmt<5) cout << ipmt << "\t" << _qfit_mpv[ipmt] << endl;
    }
    delete cdbttree;
  }
  else if ( dbase_file.EndsWith(".calib") )       // read from database
  {
    ifstream infile( dbase_location.c_str() );
    if ( !infile.is_open() )
    {
      cout << PHWHERE << "unable to open " << dbase_location << endl;
      _status = -3;
      return _status;
    }

    int pmt = -1;
    while ( infile >> pmt )
    {
      infile >> _qfit_integ[pmt] >> _qfit_mpv[pmt] >> _qfit_sigma[pmt]
        >> _qfit_integerr[pmt] >> _qfit_mpverr[pmt] >> _qfit_sigmaerr[pmt]
        >> _qfit_chi2ndf[pmt];
      if ( pmt<5 || pmt>=BBC_N_PMT-5 )
      {
        cout << pmt << "\t" <<  _qfit_integ[pmt] << "\t" <<  _qfit_mpv[pmt] << "\t" <<  _qfit_sigma[pmt]
          << "\t" <<  _qfit_integerr[pmt] << "\t" <<  _qfit_mpverr[pmt] << "\t" <<  _qfit_sigmaerr[pmt]
          << "\t" <<  _qfit_chi2ndf[pmt] << endl;
      }
    }
  }
  else
  {
    cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << endl;
    _status = -1;
    return _status;
  }

  return 1;
}

int BbcCalib::Download_TQT0(const string& dbase_location)
{
  // Reset All Values
  for (int ipmt=0; ipmt<BBC_N_PMT; ipmt++)
  {
    _tqfit_t0mean[ipmt] = NAN;
    _tqfit_t0meanerr[ipmt] = NAN;
    _tqfit_t0sigma[ipmt] = NAN;
    _tqfit_t0sigmaerr[ipmt] = NAN;
  }

  cout << "Opening " << dbase_location << endl;
  TString dbase_file = dbase_location;
  if ( dbase_file.EndsWith(".root") )       // read from database
  {
    CDBTTree *cdbttree = new CDBTTree( dbase_location );
    cdbttree->LoadCalibrations();

    for (int ipmt=0; ipmt<BBC_N_PMT; ipmt++)
    {
      _tqfit_t0mean[ipmt] = cdbttree->GetFloatValue(ipmt,"tqfit_t0mean");
      _tqfit_t0meanerr[ipmt] = cdbttree->GetFloatValue(ipmt,"tqfit_t0meanerr");
      _tqfit_t0sigma[ipmt] = cdbttree->GetFloatValue(ipmt,"tqfit_t0sigma");
      _tqfit_t0sigmaerr[ipmt] = cdbttree->GetFloatValue(ipmt,"tqfit_t0sigmaerr");
      if (ipmt<5 || ipmt>=BBC_N_PMT-5) cout << ipmt << "\t" << _tqfit_t0mean[ipmt] << endl;
    }
    delete cdbttree;
  }
  else if ( dbase_file.EndsWith(".calib") )       // read from database
  {
    ifstream infile( dbase_location.c_str() );
    if ( !infile.is_open() )
    {
      cout << PHWHERE << "unable to open " << dbase_location << endl;
      _status = -3;
      return _status;
    }

    int pmt = -1;
    while ( infile >> pmt )
    {
      infile >> _tqfit_t0mean[pmt] >> _tqfit_t0meanerr[pmt]
        >> _tqfit_t0sigma[pmt] >> _tqfit_t0sigmaerr[pmt];

      if ( pmt<5 || pmt>=BBC_N_PMT-5 )
      {
        cout << pmt << "\t" <<  _tqfit_t0mean[pmt] << "\t" <<  _tqfit_t0meanerr[pmt]
          << "\t" <<  _tqfit_t0sigma[pmt] << "\t" <<  _tqfit_t0sigmaerr[pmt] << endl;
      }
    }
    infile.close();
  }
  else
  {
    cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << endl;
    _status = -1;
    return _status;
  }

  return 1;
}

int BbcCalib::Download_SampMax(const string& dbase_location)
{
  // Reset All Values
  for (int ifeech=0; ifeech<BBC_N_FEECH; ifeech++)
  {
    _sampmax[ifeech] = -1; 
  }

  cout << "Opening " << dbase_location << endl;
  TString dbase_file = dbase_location;
  if ( dbase_file.EndsWith(".root") )       // read from database
  {
    CDBTTree *cdbttree = new CDBTTree( dbase_location );
    cdbttree->LoadCalibrations();

    for (int ifeech=0; ifeech<BBC_N_FEECH; ifeech++)
    {
      _sampmax[ifeech] = cdbttree->GetIntValue(ifeech,"sampmax");
      if (ifeech%8==0 ) cout << ifeech << "\t" << _sampmax[ifeech] << endl;
    }
    delete cdbttree;
  }
  else if ( dbase_file.EndsWith(".calib") ) // read from text file
  {
    ifstream infile( dbase_location.c_str() );
    if ( !infile.is_open() )
    {
      cout << PHWHERE << "unable to open " << dbase_location << endl;
      _status = -3;
      return _status;
    }

    int feech = -1;
    while ( infile >> feech )
    {
      infile >> _sampmax[feech];

      if ( feech<5 || feech>=BBC_N_FEECH-5 )
      {
        cout << feech << "\t" << _sampmax[feech] << endl;
      }
    }
    infile.close();
  }
  else
  {
    cout << PHWHERE << ", ERROR, unknown file type, " << dbase_location << endl;
    _status = -1;
    return _status;  // file not found
  }

  return 1;
}

int BbcCalib::StoreInDatabase()
{
  return 1;
}

void BbcCalib::Reset()
{
  // Set all initial values
  for (int ipmt=0; ipmt<BBC_N_PMT; ipmt++)
  {
    _qfit_integ[ ipmt ] = NAN;
    _qfit_mpv[ ipmt ] = NAN;
    _qfit_sigma[ ipmt ] = NAN;
    _qfit_integerr[ ipmt ] = NAN;
    _qfit_mpverr[ ipmt ] = NAN;
    _qfit_sigmaerr[ ipmt ] = NAN;
    _qfit_chi2ndf[ ipmt ] = NAN;

    _tqfit_t0mean[ ipmt ] = NAN;
    _tqfit_t0meanerr[ ipmt ] = NAN;
    _tqfit_t0sigma[ ipmt ] = NAN;
    _tqfit_t0sigmaerr[ ipmt ] = NAN;
  }

  for (int ifeech=0; ifeech<BBC_N_FEECH; ifeech++)
  {
    _sampmax[ ifeech ] = -1;
  }
}

/*
   void BbcCalib::Dump_to_file(const std::string &what)
   {
// make timerange string
string timerange = ".";
timerange.append( getStartTime()->formatTimeString() );
timerange.append("-");
timerange.append( getEndTime()->formatTimeString() );

if ( what == "ALL" || what == "GAINS" )
{
string full_outfname = "BbcCal.gains";
full_outfname.append( timerange );
ofstream outfile(full_outfname.c_str());

for (int ifeech=0; ifeech<MAXCH; ifeech++)
{
outfile << ifeech << "\t"
<< adc_gain[ifeech].getPeakChannel() << "\t"
<< adc_gain[ifeech].getDeviation() << "\t"
<< adc_gain[ifeech].getStatus() << endl;
}

outfile.close();
}

}
*/

/*
   void BbcCalib::Print(Option_t *option) const
   {
   std::string what(option);
   if ( what == "ALL" || what == "GAINS" )
   {
   cout << "GAINS " << endl;
   for (int ifeech=0; ifeech<MAXCH; ifeech++)
   {
   cout << ifeech << "\t" << adc_gain[ifeech].getPeakChannel() << endl;
   }
   }
   }
   */

