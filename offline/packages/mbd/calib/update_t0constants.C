// 
// update the t0 constants to correct the z-vertex
//
#include <mbd/MbdCalib.h>

#if defined(__CLING__)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libmbd_io.so)
#endif


void update_t0constants(const char *calib_fname, const float dz)
{
  TString cal_fname = calib_fname;

  MbdCalib *mcal = new MbdCalib();
  
  if ( cal_fname.Contains("tt_t0") )
  {
    mcal->Download_TTT0( cal_fname.Data() );
    mcal->Update_TTT0( dz );
    cal_fname.ReplaceAll( "mbd_tt_t0", "new_mbd_tt_t0" );
    mcal->Write_TTT0( cal_fname.Data() );
  }
  else if ( cal_fname.Contains("tq_t0") )
  {
    mcal->Download_TQT0( cal_fname.Data() );
    mcal->Update_TQT0(dz);
    cal_fname.ReplaceAll( "mbd_tq_t0", "new_mbd_tq_t0" );
    mcal->Write_TQT0( cal_fname.Data() );
  }
  else
  {
    cerr << "File " << calib_fname << " is not a valid t0 calibration file" << endl;
  }

}
