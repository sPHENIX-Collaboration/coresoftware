#include "TrkrNtuplizer.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase_historic/TrackSeedContainer_v1.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1RawHit.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>

#include <micromegas/MicromegasDefs.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackermillepedealignment/HelicalFitter.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

/*#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>
#include <odbc++/types.h>
*/
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>  // for shared_ptr
#include <set>     // for _Rb_tree_cons...
#include <utility>
#include <vector>

enum n_event  // NOLINT(readability-enum-initial-value, performance-enum-size)
{
  evnev,
  evnseed,
  evnrun,
  evnseg,
  evnjob,
  evsize = evnjob + 1,
};

enum n_info  // NOLINT(readability-enum-initial-value, performance-enum-size)
{
  infonocc11,
  infonocc116,
  infonocc21,
  infonocc216,
  infonocc31,
  infonocc316,
  infonrawzdc,
  infonlivezdc,
  infonscaledzdc,
  infonrawmbd,
  infonlivembd,
  infonscaledmbd,
  infonrawmbdv10,
  infonlivembdv10,
  infonscaledmbdv10,
  infonrawzdc1,
  infonlivezdc1,
  infonscaledzdc1,
  infonrawmbd1,
  infonlivembd1,
  infonscaledmbd1,
  infonrawmbdv101,
  infonlivembdv101,
  infonscaledmbdv101,
  inforzdc,
  informbd,
  informbdv10,
  infonbco1,
  infonbco,
  infonbcotr,
  infonbcotr1,
  infontrk,
  infonntpcseed,
  infonnsiseed,
  infonhitmvtx,
  infonhitintt,
  infonhittpot,
  infonhittpcall,
  infonhittpcin,
  infonhittpcmid,
  infonhittpcout,
  infonclusall,
  infonclustpc,
  infonclustpcpos,
  infonclustpcneg,
  infonclusintt,
  infonclusmvtx,
  infonclustpot,
  infosize = infonclustpot + 1
};

enum n_vertex  // NOLINT(readability-enum-initial-value, performance-enum-size)
{
  vtxnvertexID,
  vtxnvx,
  vtxnvy,
  vtxnvz,
  vtxnntracks,
  vtxnchi2,
  vtxnndof,
  vtxsize = vtxnndof + 1
};

enum n_hit  // NOLINT(readability-enum-initial-value, performance-enum-size)
{
  nhitID,
  nhite,
  nhitadc,
  nhitlayer,
  nhitphielem,
  nhitzelem,
  nhitcellID,
  nhitecell,
  nhitphibin,
  nhittbin,
  nhitphi,
  nhitr,
  nhitx,
  nhity,
  nhitz,
  hitsize = nhitz + 1
};

enum n_seed  // NOLINT(readability-enum-initial-value, performance-enum-size)
{
  nseedtrackID,
  nseedniter,
  nseedpt,
  nseedptot,
  nseedeta,
  nseedphi,
  nseedsyxint,
  nseedsrzint,
  nseedsxyslope,
  nseedsrzslope,
  nseedX0,
  nseedY0,
  nseedZ0,
  nseedR0,
  nseedcharge,
  nseeddedx,
  nseedpidedx,
  nseedkdedx,
  nseedprdedx,
  nseedn1pix,
  nseednsil,
  nseedntpc,
  nseednhits,
  seedsize = nseednhits + 1
};

enum n_residual  // NOLINT(readability-enum-initial-value, performance-enum-size)
{
  nresalpha,
  nresbeta,
  nresphio,
  nresphi,
  nresz,
  ressize = nresz + 1
};

enum n_track  // NOLINT(readability-enum-initial-value, performance-enum-size)
{
  ntrktrackID,
  ntrkcrossing,
  ntrkpx,
  ntrkpy,
  ntrkpz,
  ntrkpt,
  ntrketa,
  ntrkphi,
  ntrkdeltapt,
  ntrkdeltaeta,
  ntrkdeltaphi,
  ntrkcharge,
  ntrkquality,
  ntrkchisq,
  ntrkndf,
  ntrknhits,
  ntrknmaps,
  ntrknintt,
  ntrkntpc,
  ntrknmms,
  ntrkntpc1,
  ntrkntpc11,
  ntrkntpc2,
  ntrkntpc3,
  ntrkndedx,
  ntrknpidedx,
  ntrknkdedx,
  ntrknprdedx,
  ntrkvertexID,
  ntrkvx,
  ntrkvy,
  ntrkvz,
  ntrkdca2d,
  ntrkdca2dsigma,
  ntrkdca3dxy,
  ntrkdca3dxysigma,
  ntrkdca3dz,
  ntrkdca3dzsigma,
  ntrkpcax,
  ntrkpcay,
  ntrkpcaz,
  ntrkhlxpt,
  ntrkhlxeta,
  ntrkhlxphi,
  ntrkhlxX0,
  ntrkhlxY0,
  ntrkhlxZ0,
  ntrkhlxcharge,
  trksize = ntrkhlxcharge + 1
};

enum n_cluster  // NOLINT(readability-enum-initial-value, performance-enum-size)
{
  nclulocx,
  nclulocy,
  nclux,
  ncluy,
  ncluz,
  nclur,
  ncluphi,
  nclueta,
  nclutheta,
  ncluphibin,
  nclutbin,
  nclufee,
  ncluchan,
  nclusampa,
  ncluex,
  ncluey,
  ncluez,
  ncluephi,
  nclupez,
  nclupephi,
  nclue,
  ncluadc,
  nclumaxadc,
  ncluthick,
  ncluafac,
  nclubfac,
  ncludcal,
  nclulayer,
  ncluphielem,
  ncluzelem,
  nclusize,
  ncluphisize,
  ncluzsize,
  nclupedge,
  ncluredge,
  ncluovlp,
  nclutrackID,
  ncluniter,
  clusize = ncluniter + 1
};

TrkrNtuplizer::TrkrNtuplizer(const std::string& /*name*/, const std::string& filename, const std::string& trackmapname,
                             unsigned int nlayers_maps,
                             unsigned int nlayers_intt,
                             unsigned int nlayers_tpc,
                             unsigned int nlayers_mms)
  : SubsysReco("TrkrNtuplizer")
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _nlayers_mms(nlayers_mms)
  , _filename(filename)
  , _trackmapname(trackmapname)
{
}

TrkrNtuplizer::~TrkrNtuplizer()
{
  delete _timer;
}

int TrkrNtuplizer::Init(PHCompositeNode* /*unused*/)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");
  _tfile->SetCompressionLevel(7);

  std::string str_vertex = {"vertexID:vx:vy:vz:ntracks:chi2:ndof"};
  std::string str_event = {"event:seed:run:seg:job"};
  std::string str_hit = {"hitID:e:adc:layer:phielem:zelem:cellID:ecell:phibin:tbin:phi:r:x:y:z"};
  std::string str_cluster = {"locx:locy:x:y:z:r:phi:eta:theta:phibin:tbin:fee:chan:sampa:ex:ey:ez:ephi:pez:pephi:e:adc:maxadc:thick:afac:bfac:dcal:layer:phielem:zelem:size:phisize:zsize:pedge:redge:ovlp:trackID:niter"};
  std::string str_seed = {"seedID:siter:spt:sptot:seta:sphi:syxint:srzint:sxyslope:srzslope:sX0:sY0:sdZ0:sR0:scharge:sdedx:spidedx:skdedx:sprdedx:sn1pix:snsil:sntpc:snhits"};
  std::string str_residual = {"alpha:beta:resphio:resphi:resz"};
  std::string str_track = {"trackID:crossing:px:py:pz:pt:eta:phi:deltapt:deltaeta:deltaphi:charge:quality:chisq:ndf:nhits:nmaps:nintt:ntpc:nmms:ntpc1:ntpc11:ntpc2:ntpc3:dedx:pidedx:kdedx:prdedx:nlmaps:nlintt:nltpc:nlmms:layers:vertexID:vx:vy:vz:dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:hlxpt:hlxeta:hlxphi:hlxX0:hlxY0:hlxZ0:hlxcharge"};
  std::string str_info = {"occ11:occ116:occ21:occ216:occ31:occ316:rawzdc:livezdc:scaledzdc:rawmbd:livembd:scaledmbd:rawmbdv10:livembdv10:scaledmbdv10:rawzdc1:livezdc1:scaledzdc1:rawmbd1:livembd1:scaledmbd1:rawmbdv101:livembdv101:scaledmbdv101:rzdc:rmbd:rmbdv10:bco1:bco:bcotr:bcotr1:ntrk:ntpcseed:nsiseed:nhitmvtx:nhitintt:nhittpot:nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclustpcpos:nclustpcneg:nclusintt:nclusmaps:nclusmms"};
  if (_do_info_eval)
  {
    std::string ntp_varlist_info = str_event + ":" + str_info;
    _ntp_info = new TNtuple("ntp_info", "event info", ntp_varlist_info.c_str());
  }

  if (_do_vertex_eval)
  {
    std::string ntp_varlist_vtx = str_event + ":" + str_vertex + ":" + str_info;
    _ntp_vertex = new TNtuple("ntp_vertex", "vertex => max truth", ntp_varlist_vtx.c_str());
  }

  if (_do_hit_eval)
  {
    std::string ntp_varlist_ev = str_event + ":" + str_hit + ":" + str_info;
    _ntp_hit = new TNtuple("ntp_hit", "svtxhit => max truth", ntp_varlist_ev.c_str());
  }

  if (_do_cluster_eval)
  {
    std::string ntp_varlist_clu = str_event + ":" + str_cluster + ":" + str_info;
    _ntp_cluster = new TNtuple("ntp_cluster", "svtxcluster => max truth", ntp_varlist_clu.c_str());
  }
  if (_do_clus_trk_eval)
  {
    std::string ntp_varlist_clut = str_event + ":" + str_cluster + ":" + str_residual + ":" + str_seed + ":" + str_info;
    _ntp_clus_trk = new TNtuple("ntp_clus_trk", "cluster on track", ntp_varlist_clut.c_str());
  }

  if (_do_track_eval)
  {
    std::string ntp_varlist_trk = str_event + ":" + str_track + ":" + str_info;
    _ntp_track = new TNtuple("ntp_track", "svtxtrack => max truth", ntp_varlist_trk.c_str());
  }

  if (_do_tpcseed_eval)
  {
    std::string ntp_varlist_tsee = str_event + ":" + str_seed + ":" + str_info;
    _ntp_tpcseed = new TNtuple("ntp_tpcseed", "seeds from truth", ntp_varlist_tsee.c_str());
  }
  if (_do_siseed_eval)
  {
    std::string ntp_varlist_ssee = str_event + ":" + str_seed + ":" + str_info;
    _ntp_siseed = new TNtuple("ntp_siseed", "seeds from truth", ntp_varlist_ssee.c_str());
  }

  std::string dedx_fitparams = CDBInterface::instance()->getUrl("TPC_DEDX_FITPARAM");
  TFile* filefit = TFile::Open(dedx_fitparams.c_str());

  if (!filefit->IsOpen())
  {
    std::cout << "Error opening filefit!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  filefit->GetObject("f_piband", f_pion_plus);
  filefit->GetObject("f_Kband", f_kaon_plus);
  filefit->GetObject("f_pband", f_proton_plus);
  filefit->GetObject("f_piminus_band", f_pion_minus);
  filefit->GetObject("f_Kminus_band", f_kaon_minus);
  filefit->GetObject("f_pbar_band", f_proton_minus);

  _timer = new PHTimer("_eval_timer");
  _timer->stop();
  /**/

  dedxcorr[0][0][0] = 0.54;
  dedxcorr[0][0][1] = 0.90;
  dedxcorr[0][0][2] = 1.11;
  dedxcorr[0][1][0] = 0.72;
  dedxcorr[0][1][1] = 0.98;
  dedxcorr[0][1][2] = 0.99;
  dedxcorr[0][2][0] = 0.83;
  dedxcorr[0][2][1] = 0.81;
  dedxcorr[0][2][2] = 1.03;
  dedxcorr[0][3][0] = 0.82;
  dedxcorr[0][3][1] = 0.52;
  dedxcorr[0][3][2] = 1.02;
  dedxcorr[0][4][0] = 0.85;
  dedxcorr[0][4][1] = 0.78;
  dedxcorr[0][4][2] = 0.82;
  dedxcorr[0][5][0] = 0.85;
  dedxcorr[0][5][1] = 0.70;
  dedxcorr[0][5][2] = 0.92;
  dedxcorr[0][6][0] = 0.78;
  dedxcorr[0][6][1] = 0.94;
  dedxcorr[0][6][2] = 1.03;
  dedxcorr[0][7][0] = 0.64;
  dedxcorr[0][7][1] = 0.94;
  dedxcorr[0][7][2] = 1.19;
  dedxcorr[0][8][0] = 0.95;
  dedxcorr[0][8][1] = 0.47;
  dedxcorr[0][8][2] = 1.12;
  dedxcorr[0][9][0] = 0.74;
  dedxcorr[0][9][1] = 1.07;
  dedxcorr[0][9][2] = 0.85;
  dedxcorr[0][10][0] = 0.86;
  dedxcorr[0][10][1] = 0.74;
  dedxcorr[0][10][2] = 1.08;
  dedxcorr[0][11][0] = 0.60;
  dedxcorr[0][11][1] = 1.03;
  dedxcorr[0][11][2] = 0.89;
  dedxcorr[1][0][0] = 0.74;
  dedxcorr[1][0][1] = 0.76;
  dedxcorr[1][0][2] = 0.89;
  dedxcorr[1][1][0] = 0.61;
  dedxcorr[1][1][1] = 0.87;
  dedxcorr[1][1][2] = 1.10;
  dedxcorr[1][2][0] = 0.73;
  dedxcorr[1][2][1] = 0.86;
  dedxcorr[1][2][2] = 1.70;
  dedxcorr[1][3][0] = 0.40;
  dedxcorr[1][3][1] = 0.76;
  dedxcorr[1][3][2] = 0.99;
  dedxcorr[1][4][0] = 0.82;
  dedxcorr[1][4][1] = 0.71;
  dedxcorr[1][4][2] = 0.86;
  dedxcorr[1][5][0] = 0.65;
  dedxcorr[1][5][1] = 1.00;
  dedxcorr[1][5][2] = 0.79;
  dedxcorr[1][6][0] = 0.75;
  dedxcorr[1][6][1] = 0.84;
  dedxcorr[1][6][2] = 0.40;
  dedxcorr[1][7][0] = 0.95;
  dedxcorr[1][7][1] = 0.81;
  dedxcorr[1][7][2] = 0.85;
  dedxcorr[1][8][0] = 0.79;
  dedxcorr[1][8][1] = 0.67;
  dedxcorr[1][8][2] = 1.17;
  dedxcorr[1][9][0] = 0.73;
  dedxcorr[1][9][1] = 0.94;
  dedxcorr[1][9][2] = 0.95;
  dedxcorr[1][10][0] = 0.40;
  dedxcorr[1][10][1] = 0.83;
  dedxcorr[1][10][2] = 0.87;
  dedxcorr[1][11][0] = 0.80;
  dedxcorr[1][11][1] = 0.63;
  dedxcorr[1][11][2] = 0.68;
  dedxcorr[0][0][1] *= 0.88;
  dedxcorr[0][0][2] *= 0.92;
  dedxcorr[0][1][0] *= 0.90;
  dedxcorr[0][1][1] *= 0.92;
  dedxcorr[0][1][2] *= 0.94;
  dedxcorr[0][2][0] *= 0.80;
  dedxcorr[0][2][1] *= 0.81;
  dedxcorr[0][2][2] *= 0.88;
  dedxcorr[0][3][0] *= 0.76;
  dedxcorr[0][3][1] *= 0.76;
  dedxcorr[0][3][2] *= 0.77;
  dedxcorr[0][4][0] *= 0.87;
  dedxcorr[0][4][1] *= 0.88;
  dedxcorr[0][4][2] *= 0.85;
  dedxcorr[0][5][0] *= 0.86;
  dedxcorr[0][5][1] *= 0.87;
  dedxcorr[0][5][2] *= 0.92;
  dedxcorr[0][6][0] *= 0.89;
  dedxcorr[0][6][1] *= 0.91;
  dedxcorr[0][6][2] *= 0.91;
  dedxcorr[0][7][0] *= 0.90;
  dedxcorr[0][7][1] *= 0.90;
  dedxcorr[0][7][2] *= 0.94;
  dedxcorr[0][8][0] *= 0.87;
  dedxcorr[0][8][1] *= 0.90;
  dedxcorr[0][8][2] *= 0.91;
  dedxcorr[0][9][0] *= 0.92;
  dedxcorr[0][9][1] *= 0.92;
  dedxcorr[0][9][2] *= 0.94;
  dedxcorr[0][10][0] *= 0.91;
  dedxcorr[0][10][1] *= 0.92;
  dedxcorr[0][10][2] *= 0.93;
  dedxcorr[0][11][0] *= 0.87;
  dedxcorr[0][11][1] *= 0.88;
  dedxcorr[0][11][2] *= 0.91;
  dedxcorr[1][0][0] *= 0.91;
  dedxcorr[1][0][1] *= 0.92;
  dedxcorr[1][0][2] *= 0.92;
  dedxcorr[1][1][0] *= 0.92;
  dedxcorr[1][1][1] *= 0.89;
  dedxcorr[1][1][2] *= 0.91;
  dedxcorr[1][2][0] *= 0.84;
  dedxcorr[1][2][1] *= 0.87;
  dedxcorr[1][2][2] *= 1.05;
  dedxcorr[1][3][0] *= 0.83;
  dedxcorr[1][3][1] *= 0.89;
  dedxcorr[1][3][2] *= 0.90;
  dedxcorr[1][4][0] *= 0.90;
  dedxcorr[1][4][1] *= 0.91;
  dedxcorr[1][4][2] *= 0.91;
  dedxcorr[1][5][0] *= 0.92;
  dedxcorr[1][5][1] *= 0.93;
  dedxcorr[1][5][2] *= 0.93;
  dedxcorr[1][6][0] *= 0.87;
  dedxcorr[1][6][1] *= 0.85;
  dedxcorr[1][6][2] *= 0.99;
  dedxcorr[1][7][0] *= 0.92;
  dedxcorr[1][7][1] *= 0.89;
  dedxcorr[1][7][2] *= 0.89;
  dedxcorr[1][8][0] *= 0.88;
  dedxcorr[1][8][1] *= 0.89;
  dedxcorr[1][8][2] *= 0.95;
  dedxcorr[1][9][0] *= 0.89;
  dedxcorr[1][9][1] *= 0.90;
  dedxcorr[1][9][2] *= 0.91;
  dedxcorr[1][10][0] *= 0.86;
  dedxcorr[1][10][1] *= 0.92;
  dedxcorr[1][10][2] *= 0.95;
  dedxcorr[1][11][0] *= 0.85;
  dedxcorr[1][11][1] *= 0.85;
  dedxcorr[1][11][2] *= 0.88;

  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrNtuplizer::InitRun(PHCompositeNode* topNode)
{
  auto* geom =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  AdcClockPeriod = geom->GetFirstLayerCellGeom()->get_zstep();

  // Create Fee Map
  auto* geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  {
    if (!geom_container)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  m_cdb = CDBInterface::instance();
  std::string calibdir = m_cdb->getUrl("TPC_FEE_CHANNEL_MAP");

  if (calibdir[0] == '/')
  {
    // use generic CDBTree to load
    m_cdbttree = new CDBTTree(calibdir);
    m_cdbttree->LoadCalibrations();
  }
  else
  {
    std::cout << "TrkrNtuplizer::::InitRun No calibration file found" << std::endl;
    exit(1);
  }
  //  for(unsigned int side = 0; side <=1;side++){

  for (unsigned int sector = 0; sector < 24; sector++)
  {
    int side = 1;
    if (sector > 11)
    {
      side = 0;
    }
    for (unsigned int fee = 0; fee <= 25; fee++)
    {
      for (unsigned int channel = 0; channel <= 255; channel++)
      {
        int feeM = FEE_map[fee];
        if (FEE_R[fee] == 2)
        {
          feeM += 6;
        }
        if (FEE_R[fee] == 3)
        {
          feeM += 14;
        }
        unsigned int key = (256 * (feeM)) + channel;
        std::string varname = "layer";
        unsigned int layer = m_cdbttree->GetIntValue(key, varname);
        // antenna pads will be in 0 layer
        if (layer <= 0 or layer>55)
        {
          continue;
        }
        varname = "phi";  // + std::to_string(key);
        double phi = -1 * pow(-1, side) * m_cdbttree->GetDoubleValue(key, varname) + (sector % 12) * M_PI / 6;
        PHG4TpcCylinderGeom* layergeom = geom_container->GetLayerCellGeom(layer);
        unsigned int phibin = layergeom->get_phibin(phi);
        // get global coords
        double radius = layergeom->get_radius();  // returns center of layer
        double chanphi = layergeom->get_phi(phibin);
        float chanx = radius * cos(chanphi);
        float chany = radius * sin(chanphi);

        TrkrDefs::cluskey dummykey = TpcDefs::genClusKey(layer, (mc_sectors[sector % 12]), side, phibin);
        fee_info nufeeinfo;
        nufeeinfo.fee = fee;
        nufeeinfo.channel = channel;
        nufeeinfo.sampa = channel % 32;
        fee_map.insert(std::make_pair(dummykey, nufeeinfo));
        if (Verbosity() > 0)
        {
          std::cout << "side sec lay phibin " << side
                    << " | hw " << sector
                    << " | mc " << (mc_sectors[sector % 12])
                    << " | " << layer
                    << " | " << phibin
                    << " => fee, channel "
                    << fee << " | " << channel
                    << " x | y " << chanx << " | " << chany
                    << std::endl;
        }
        // have side, sector, layer, phibin -> fee, channel
        //	  hit_set_key = TpcDefs::genHitSetKey(layer, (mc_sectors[sector % 12]), side);
        //   hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrNtuplizer::process_event(PHCompositeNode* topNode)
{
  /*
  cout << "HELLO" << endl;
  //get ZDC rate
  odbc::Connection *dbConnection = odbc::DriverManager::getConnection("daq", "", "");
  std::string sql = "SELECT * FROM gl1_scalers WHERE runnumber = " + std::to_string(m_runnumber) + ";";
  odbc::Statement *stmt = dbConnection->createStatement();
  odbc::ResultSet *resultSet = stmt->executeQuery(sql);
  std::array<std::array<uint64_t, 3>, 64> scalers{};  // initialize to zero
  if (!resultSet)
  {
    std::cerr << "No db found for run number " << m_runnumber << ". Cannot get ZDC rate so aborting run" << std
::endl;
    delete resultSet;
    delete stmt;
    delete dbConnection;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  while (resultSet->next())
  {
    int index = resultSet->getInt("index");
    // Iterate over the columns and fill the TriggerRunInfo object
    scalers[index][0] = resultSet->getLong("scaled");
    scalers[index][1] = resultSet->getLong("live");
    scalers[index][2] = resultSet->getLong("raw");
  }

  delete resultSet;
  delete stmt;
  delete dbConnection;
  m_ZDC_coincidence = (1.0*scalers[3][2]/scalers[0][2])/(106e-9);

  MbdPmtContainer* bbcpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  m_mbd_rate = 0;
  if (bbcpmts)
  {
    int nPMTs = bbcpmts->get_npmt();
    for (int i = 0; i < nPMTs; i++)
    {
      m_mbd_rate += bbcpmts->get_pmt(i)->get_q();
    }
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << "TrackResiduals::process_event: Could not find MbdPmtContainer," << std::endl;
    }
  }
  */
  auto* gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
  if (gl1)
  {
    m_bco = gl1->get_bco();
    auto lbshift = m_bco << 24U;
    m_bcotr = lbshift >> 24U;
  }
  else
  {
    Gl1Packet* gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
    if (!gl1PacketInfo)
    {
      m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
      m_bcotr = std::numeric_limits<uint64_t>::quiet_NaN();
    }
    m_firedTriggers.clear();

    if (gl1PacketInfo)
    {
      /*      uint64_t raw_scalers[64] = {0};
      uint64_t live_scalers[64] = {0};
      uint64_t scaled_scalers[64] = {0};
      */
      m_gl1BunchCrossing = gl1PacketInfo->getBunchNumber();
      uint64_t triggervec = gl1PacketInfo->getTriggerVector();
      m_bco = gl1PacketInfo->getBCO();
      auto lbshift = m_bco << 24U;
      m_bcotr = lbshift >> 24U;
      for (int i = 0; i < 64; i++)
      {
        /*
        raw_scalers[i] =  gl1PacketInfo->lValue(0, i);
        live_scalers[i] =   gl1PacketInfo->lValue(1, i);
        scaled_scalers[i] =  gl1PacketInfo->lValue(2, i);
        */
        bool trig_decision = ((triggervec & 0x1U) == 0x1U);
        if (trig_decision)
        {
          m_firedTriggers.push_back(i);
        }
        triggervec = (triggervec >> 1U) & 0xffffffffU;
      }
      for (int i1 = 0; i1 <= 12; i1++)
      {
        for (int i2 = 0; i2 < 3; i2++)
        {
          //	  if (i1>=3&&i1<=9) continue;
          // if (i1==11) continue;
          // if (i1==10&&i2==2) continue;
          /*
          std::cout << " index 1: " << i1
                    << " index 2: " << i2
                    << " scaler: " << gl1PacketInfo->lValue(i1, i2)
                    << std::endl;
          */
        }
      }
      m_rawzdc = gl1PacketInfo->lValue(2, 0);         // raw_scalers[1];
      m_livezdc = gl1PacketInfo->lValue(2, 1);        // live_scalers[1];
      m_scaledzdc = gl1PacketInfo->lValue(2, 2);      // scaled_scalers[1];
      m_rawmbd = gl1PacketInfo->lValue(10, 0);        // raw_scalers[10];
      m_livembd = gl1PacketInfo->lValue(10, 1);       // live_scalers[10];
      m_scaledmbd = gl1PacketInfo->lValue(10, 2);     // scaled_scalers[10];
      m_rawmbdv10 = gl1PacketInfo->lValue(12, 0);     // raw_scalers[12];
      m_livembdv10 = gl1PacketInfo->lValue(12, 1);    // live_scalers[12];
      m_scaledmbdv10 = gl1PacketInfo->lValue(12, 2);  // scaled_scalers[12];

      m_rawzdc_hist.push(m_rawzdc);
      m_rawmbd_hist.push(m_rawmbd);
      m_rawmbdv10_hist.push(m_rawmbdv10);
      m_bco_hist.push(m_bco);

      /*    MBD NS >= 2 is trigger bit 10
            MBD NS >=2 w/ vtx < 10 cm is trigger bit 12
            ZDC is trigger bit 1
      */
    }
  }

  if ((Verbosity() > 1) && (_ievent % 100 == 0))
  {
    std::cout << "TrkrNtuplizer::process_event - Event = " << _ievent << std::endl;
  }

  recoConsts* rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
  {
    _iseed = rc->get_IntFlag("RANDOMSEED");
    m_fSeed = _iseed;
  }
  else
  {
    _iseed = 0;
    m_fSeed = std::numeric_limits<float>::quiet_NaN();
  }
  if (_trackmap == nullptr)
  {
    _trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname);
  }

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if (!_cluster_map)
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }

  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tgeometry)
  {
    std::cout << "No Acts geometry on node tree. Can't  continue."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  {
    if (!_geom_container)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << "TrkrNtuplizer::process_event - Seed = " << _iseed << std::endl;
  }
  //-----------------------------------
  // print what is coming into the code
  //-----------------------------------

  printInputInfo(topNode);

  //---------------------------
  // fill the Evaluator NTuples
  //---------------------------

  fillOutputNtuples(topNode);

  //--------------------------------------------------
  // Print out the ancestry information for this event
  //--------------------------------------------------

  // printOutputInfo(topNode);
  /*
  ++_ievent;
  if(m_rawzdc_hist.size()==50){
    m_rawzdc_hist.pop();
    m_rawmbd_hist.pop();
    m_rawmbdv10_hist.pop();
    m_bco_hist.pop();
  }
  m_rawzdclast = m_rawzdc;
  m_rawmbdlast = m_rawmbd;
  m_rawmbdv10last = m_rawmbdv10;
  m_bcolast =m_bco;
  */
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrNtuplizer::End(PHCompositeNode* /*topNode*/)
{
  _tfile->cd();

  if (_ntp_info)
  {
    _ntp_info->Write();
  }
  if (_ntp_vertex)
  {
    _ntp_vertex->Write();
  }
  if (_ntp_hit)
  {
    _ntp_hit->Write();
  }
  if (_ntp_cluster)
  {
    _ntp_cluster->Write();
  }
  if (_ntp_clus_trk)
  {
    _ntp_clus_trk->Write();
  }
  if (_ntp_track)
  {
    _ntp_track->Write();
  }
  if (_ntp_tpcseed)
  {
    _ntp_tpcseed->Write();
  }
  if (_ntp_siseed)
  {
    _ntp_siseed->Write();
  }

  _tfile->Close();

  delete _tfile;

  if (Verbosity() > 1)
  {
    std::cout << "========================= TrkrNtuplizer::End() ============================" << std::endl;
    std::cout << " " << _ievent << " events of output written to: " << _filename << std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrkrNtuplizer::printInputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "TrkrNtuplizer::printInputInfo() entered" << std::endl;
  }

  if (Verbosity() > 3)
  {
    // event information
    std::cout << std::endl;
    std::cout << PHWHERE << "   INPUT FOR EVENT " << _ievent << std::endl;

    std::cout << std::endl;

    std::cout << "---SVTXCLUSTERS-------------" << std::endl;

    if (_cluster_map != nullptr)
    {
      unsigned int icluster = 0;
      for (const auto& hitsetkey : _cluster_map->getHitSetKeys())
      {
        auto range = _cluster_map->getClusters(hitsetkey);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          TrkrDefs::cluskey cluster_key = iter->first;
          std::cout << icluster << " with key " << cluster_key << " of " << _cluster_map->size();
          std::cout << ": SvtxCluster: " << std::endl;
          iter->second->identify();
          ++icluster;
        }
      }
    }

    std::cout << "---SVXTRACKS-------------" << std::endl;

    if (_trackmap)
    {
      unsigned int itrack = 0;
      for (auto& iter : *_trackmap)
      {
        std::cout << itrack << " of " << _trackmap->size();
        SvtxTrack* track = iter.second;
        std::cout << " : SvtxTrack:" << std::endl;
        track->identify();
        std::cout << std::endl;
        ++itrack;
      }
    }

    std::cout << "---SVXVERTEXES-------------" << std::endl;
    SvtxVertexMap* vertexmap = nullptr;
    vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices

    if (vertexmap)
    {
      unsigned int ivertex = 0;
      for (SvtxVertexMap::Iter iter = vertexmap->begin();
           iter != vertexmap->end();
           ++iter)
      {
        std::cout << ivertex << " of " << vertexmap->size();
        SvtxVertex* vertex = iter->second;
        std::cout << " : SvtxVertex:" << std::endl;
        vertex->identify();
        std::cout << std::endl;
      }
    }
  }

  return;
}

void TrkrNtuplizer::printOutputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "TrkrNtuplizer::printOutputInfo() entered" << std::endl;
  }

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if (Verbosity() > 100)
  {
    // event information
    std::cout << std::endl;
    std::cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << std::endl;
    std::cout << std::endl;

    float vx = std::numeric_limits<float>::quiet_NaN();
    float vy = std::numeric_limits<float>::quiet_NaN();
    float vz = std::numeric_limits<float>::quiet_NaN();

    SvtxVertexMap* vertexmap = nullptr;

    vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices

    if (vertexmap)
    {
      if (!vertexmap->empty())
      {
        SvtxVertex* vertex = (vertexmap->begin()->second);

        vx = vertex->get_x();
        vy = vertex->get_y();
        vz = vertex->get_z();
      }
    }

    std::cout << "===Vertex Reconstruction=======================" << std::endl;
    std::cout << "vreco = (" << vx << "," << vy << "," << vz << ")" << std::endl;
    std::cout << std::endl;

    std::cout << "===Tracking Summary============================" << std::endl;

    TrkrHitSetContainer* hitsetmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

    unsigned int nclusters[100] = {0};
    unsigned int nhits[100] = {0};

    _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

    if (!_tgeometry)
    {
      std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
                << std::endl;
    }

    if (hitsetmap)
    {
      TrkrHitSetContainer::ConstRange all_hitsets = hitsetmap->getHitSets();
      for (TrkrHitSetContainer::ConstIterator hitsetiter = all_hitsets.first;
           hitsetiter != all_hitsets.second;
           ++hitsetiter)
      {
        // we have a single hitset, get the layer
        unsigned int layer = TrkrDefs::getLayer(hitsetiter->first);

        // count all hits in this hitset
        TrkrHitSet::ConstRange hitrangei = hitsetiter->second->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {
          ++nhits[layer];
        }
        auto range = _cluster_map->getClusters(hitsetiter->first);
        for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
        {
          const auto cluskey = clusIter->first;
          nclusters[TrkrDefs::getLayer(cluskey)]++;
        }
      }
    }

    for (unsigned int ilayer = 0; ilayer < _nlayers_maps + _nlayers_intt + _nlayers_tpc; ++ilayer)
    {
      PHG4TpcCylinderGeom* GeoLayer = _geom_container->GetLayerCellGeom(ilayer);

      std::cout << "layer " << ilayer
                << " => nHits = " << nhits[ilayer]
                << " => nClusters = " << nclusters[ilayer] << std::endl;
      if (ilayer >= _nlayers_maps + _nlayers_intt && ilayer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
      {
        std::cout << "layer " << ilayer
                  << " => nphi = " << GeoLayer->get_phibins()
                  << " => nz   = " << GeoLayer->get_zbins()
                  << " => ntot = " << GeoLayer->get_phibins() * GeoLayer->get_zbins()
                  << std::endl;
      }
    }

    std::cout << " => nTracks = ";
    if (_trackmap)
    {
      std::cout << _trackmap->size() << std::endl;
    }
    else
    {
      std::cout << 0 << std::endl;
    }

    std::cout << std::endl;

  }  // if Verbosity()

  return;
}

void TrkrNtuplizer::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "TrkrNtuplizer::fillOutputNtuples() entered" << std::endl;
  }

  float fx_event[n_event::evsize]{(float) _ievent, (float) _iseed, (float) m_runnumber, (float) m_segment, (float) m_job};
  float fx_info[((int) (n_info::infosize))] = {0};

  float nhit[100];
  for (float& i : nhit)
  {
    i = 0;
  }

  if (_ntp_info)
  {
    TrkrHitSetContainer* hitmap_in = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    if (hitmap_in)
    {
      TrkrHitSetContainer::ConstRange all_hitsets = hitmap_in->getHitSets();
      for (TrkrHitSetContainer::ConstIterator hitsetiter = all_hitsets.first;
           hitsetiter != all_hitsets.second;
           ++hitsetiter)
      {
        // we have a single hitset, get the layer
        unsigned int layer = TrkrDefs::getLayer(hitsetiter->first);

        TrkrHitSet::ConstRange hitrangei = hitsetiter->second->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {
          nhit[layer]++;

          if (layer < _nlayers_maps)
          {
            fx_info[n_info::infonhitmvtx]++;
          }
          if (layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt)
          {
            fx_info[n_info::infonhitintt]++;
          }
          if ((float) layer >= _nlayers_maps + _nlayers_intt &&
              (float) layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            fx_info[n_info::infonhittpcall]++;
          }
          if ((float) layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            fx_info[n_info::infonhittpot]++;
          }
          if ((float) layer == _nlayers_maps + _nlayers_intt)
          {
            fx_info[n_info::infonhittpcin]++;
          }
          if ((float) layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc - 1)
          {
            fx_info[n_info::infonhittpcout]++;
          }
          // NOLINTNEXTLINE(bugprone-integer-division)
          if ((float) layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc / 2 - 1)
          {
            fx_info[n_info::infonhittpcmid]++;
          }
        }
      }
    }

    _geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    if (!_geom_container)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return;
    }
    PHG4TpcCylinderGeom* GeoLayer;
    int layer = _nlayers_maps + _nlayers_intt;
    GeoLayer = _geom_container->GetLayerCellGeom(layer);
    int nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
    float nhits = nhit[layer];
    fx_info[n_info::infonocc11] = nhits / nbins;

    layer = _nlayers_maps + _nlayers_intt + 15;
    GeoLayer = _geom_container->GetLayerCellGeom(layer);
    nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
    nhits = nhit[layer];
    fx_info[n_info::infonocc116] = nhits / nbins;

    layer = _nlayers_maps + _nlayers_intt + 16;
    GeoLayer = _geom_container->GetLayerCellGeom(layer);
    nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
    nhits = nhit[layer];
    fx_info[n_info::infonocc21] = nhits / nbins;
    layer = _nlayers_maps + _nlayers_intt + 31;
    GeoLayer = _geom_container->GetLayerCellGeom(layer);
    nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
    nhits = nhit[layer];
    fx_info[n_info::infonocc216] = nhits / nbins;
    layer = _nlayers_maps + _nlayers_intt + 32;
    GeoLayer = _geom_container->GetLayerCellGeom(layer);
    nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
    nhits = nhit[layer];
    fx_info[n_info::infonocc31] = nhits / nbins;
    layer = _nlayers_maps + _nlayers_intt + 47;
    GeoLayer = _geom_container->GetLayerCellGeom(layer);
    nbins = GeoLayer->get_phibins() * GeoLayer->get_zbins();
    nhits = nhit[layer];
    fx_info[n_info::infonocc316] = nhits / nbins;

    if (Verbosity() > 1)
    {
      std::cout << " occ11 = " << fx_info[n_info::infonocc11]
                << " occ116 = " << fx_info[n_info::infonocc116]
                << " occ21 = " << fx_info[n_info::infonocc21]
                << " occ216 = " << fx_info[n_info::infonocc216]
                << " occ31 = " << fx_info[n_info::infonocc31]
                << " occ316 = " << fx_info[n_info::infonocc316]
                << std::endl;
    }
    /*
    std::cout << " bco = " << m_bco << " " << m_bcolast
         << " zdc = " << m_rawzdc << " " << m_rawzdclast
         << " mbd = " << m_rawmbd << " " << m_rawmbdlast
         << std::endl;

    std::cout << " zdc rate = "    << (m_rawzdc_hist.back()   -m_rawzdc_hist.front())/((m_bco_hist.back() - m_bco_hist.front())*106e-9)
         << " mbd rate = "    << (m_rawmbd_hist.back()   -m_rawmbd_hist.front())/((m_bco_hist.back() - m_bco_hist.front())*106e-9)
         << " mbdv10 rate = " << (m_rawmbdv10_hist.back()-m_rawmbdv10_hist.front())/((m_bco_hist.back() - m_bco_hist.front())*106e-9)
         << std::endl;
    */
    fx_info[n_info::inforzdc] = (m_rawzdc_hist.back() - m_rawzdc_hist.front()) / ((m_bco_hist.back() - m_bco_hist.front()) * 106e-9);
    fx_info[n_info::informbd] = (m_rawmbd_hist.back() - m_rawmbd_hist.front()) / ((m_bco_hist.back() - m_bco_hist.front()) * 106e-9);
    fx_info[n_info::informbdv10] = (m_rawmbd_hist.back() - m_rawmbd_hist.front()) / ((m_bco_hist.back() - m_bco_hist.front()) * 106e-9);
    if (_ievent == 0)
    {
      m_bco1 = m_bco;
      m_bcotr1 = m_bcotr;
      m_mbd_rate1 = m_mbd_rate;
      m_rawzdc1 = m_rawzdc;
      m_livezdc1 = m_livezdc;
      m_scaledzdc1 = m_scaledzdc;
      m_rawmbd1 = m_rawmbd;
      m_livembd1 = m_livembd;
      m_scaledmbd1 = m_scaledmbd;
      m_rawmbdv101 = m_rawmbdv10;
      m_livembdv101 = m_livembdv10;
      m_scaledmbdv101 = m_scaledmbdv10;
    }
    fx_info[n_info::infonbco] = m_bco - m_bco1;
    fx_info[n_info::infonbcotr] = m_bcotr - m_bcotr1;
    fx_info[n_info::infonbco1] = m_bco1;
    fx_info[n_info::infonbcotr1] = m_bcotr1;

    fx_info[n_info::infonrawzdc] = m_rawzdc;
    fx_info[n_info::infonlivezdc] = m_livezdc;
    fx_info[n_info::infonscaledzdc] = m_scaledzdc;
    fx_info[n_info::infonrawmbd] = m_rawmbd;
    fx_info[n_info::infonlivembd] = m_livembd;
    fx_info[n_info::infonscaledmbd] = m_scaledmbd;
    fx_info[n_info::infonrawmbdv10] = m_rawmbdv10;
    fx_info[n_info::infonlivembdv10] = m_livembdv10;
    fx_info[n_info::infonscaledmbdv10] = m_scaledmbdv10;

    fx_info[n_info::infonrawzdc1] = m_rawzdc1;
    fx_info[n_info::infonlivezdc1] = m_livezdc1;
    fx_info[n_info::infonscaledzdc1] = m_scaledzdc1;
    fx_info[n_info::infonrawmbd1] = m_rawmbd1;
    fx_info[n_info::infonlivembd1] = m_livembd1;
    fx_info[n_info::infonscaledmbd1] = m_scaledmbd1;
    fx_info[n_info::infonrawmbdv101] = m_rawmbdv101;
    fx_info[n_info::infonlivembdv101] = m_livembdv101;
    fx_info[n_info::infonscaledmbdv101] = m_scaledmbdv101;

    if (_cluster_map)
    {
      fx_info[n_info::infonclusall] = _cluster_map->size();

      for (const auto& hitsetkey : _cluster_map->getHitSetKeys())
      {
        auto range = _cluster_map->getClusters(hitsetkey);
        for (auto iter_cin = range.first; iter_cin != range.second; ++iter_cin)
        {
          TrkrDefs::cluskey cluster_key = iter_cin->first;
          unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
          unsigned int side_local = TpcDefs::getSide(cluster_key);
          if (_nlayers_maps > 0)
          {
            if (layer_local < _nlayers_maps)
            {
              fx_info[n_info::infonclusmvtx]++;
            }
          }
          if (_nlayers_intt > 0)
          {
            if ((layer_local >= _nlayers_maps) && (layer_local < (_nlayers_maps + _nlayers_intt)))
            {
              fx_info[n_info::infonclusintt]++;
            }
          }
          if (_nlayers_tpc > 0)
          {
            if (layer_local >= (_nlayers_maps + _nlayers_intt) && layer_local < (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
            {
              fx_info[n_info::infonclustpc]++;
              if (side_local == 0)
              {
                fx_info[n_info::infonclustpcneg]++;
              }
              if (side_local == 1)
              {
                fx_info[n_info::infonclustpcpos]++;
              }
            }
          }
          if (_nlayers_mms > 0)
          {
            if (layer_local >= (_nlayers_maps + _nlayers_intt + _nlayers_tpc))
            {
              fx_info[n_info::infonclustpot]++;
            }
          }
        }
      }
    }
  }
  //-----------------------
  // fill the info NTuple
  //-----------------------
  if (_ntp_info)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_info " << std::endl;
    }

    if (_trackmap)
    {
      fx_info[n_info::infontrk] = (float) _trackmap->size();
    }
    auto* siseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
    auto* tpcseedmap = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
    if (tpcseedmap)
    {
      fx_info[n_info::infonntpcseed] = tpcseedmap->size();
    }
    else
    {
      fx_info[n_info::infonntpcseed] = 0;
    }
    if (siseedmap)
    {
      fx_info[n_info::infonnsiseed] = siseedmap->size();
    }
    else
    {
      fx_info[n_info::infonnsiseed] = 0;
    }
    if (Verbosity() > 0)
    {
      std::cout << "EVENTINFO SEED: " << m_fSeed << std::endl;
      std::cout << "EVENTINFO NHIT: " << std::setprecision(9) << fx_info[n_info::infonclusall] << std::endl;
      std::cout << "EVENTINFO CLUSTPC: " << fx_info[n_info::infonclustpc] << std::endl;
      std::cout << "EVENTINFO NTRKREC: " << fx_info[n_info::infontrk] << std::endl;
    }

    float* info_data = new float[((int) (n_info::infosize)) + n_event::evsize];
    std::copy(fx_event, fx_event + n_event::evsize, info_data);
    std::copy(fx_info, fx_info + ((int) (n_info::infosize)), info_data + n_event::evsize);
    _ntp_info->Fill(info_data);
    delete[] info_data;
  }

  //-----------------------
  // fill the Vertex NTuple
  //-----------------------
  bool doit = true;
  if (_ntp_vertex && doit)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_vertex " << std::endl;
      std::cout << "start vertex time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
      _timer->restart();
    }
    float fx_vertex[n_vertex::vtxsize];
    for (float& i : fx_vertex)
    {
      i = 0;
    }

    //    SvtxVertexMap* vertexmap = nullptr;

    //    vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices

    float vx = std::numeric_limits<float>::quiet_NaN();
    float vy = std::numeric_limits<float>::quiet_NaN();
    float vz = std::numeric_limits<float>::quiet_NaN();
    float ntracks = std::numeric_limits<float>::quiet_NaN();
    fx_vertex[vtxnvx] = vx;
    fx_vertex[vtxnvy] = vy;
    fx_vertex[vtxnvz] = vz;
    fx_vertex[vtxnntracks] = ntracks;
    if (Verbosity() > 1)
    {
      std::cout << " adding vertex data " << std::endl;
    }
    float* vertex_data = new float[((int) (n_info::infosize)) + n_event::evsize + n_vertex::vtxsize];
    std::copy(fx_event, fx_event + n_event::evsize, vertex_data);
    std::copy(fx_vertex, fx_vertex + n_vertex::vtxsize, vertex_data + n_event::evsize);
    std::copy(fx_info, fx_info + ((int) (n_info::infosize)), vertex_data + n_event::evsize + n_vertex::vtxsize);
    _ntp_vertex->Fill(vertex_data);
    delete[] vertex_data;
  }
  if (Verbosity() > 1)
  {
    _timer->stop();
    std::cout << "vertex time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
  }
  //--------------------
  // fill the Hit NTuple
  //--------------------

  if (_ntp_hit)
  {
    float* hit_data = new float[((int) (n_info::infosize)) + n_event::evsize + n_hit::hitsize];
    float fx_hit_0[((int) (n_info::infosize)) + n_event::evsize + n_hit::hitsize] = {0};
    float fx_hit[((int) (n_hit::hitsize))] = {0};
    auto* m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

    if (Verbosity() >= 1)
    {
      std::cout << "Filling ntp_hit " << std::endl;
      _timer->restart();
    }
    // need things off of the DST...
    TrkrHitSetContainer* hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

    if (hitmap)
    {
      TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets();
      for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
           iter != all_hitsets.second;
           ++iter)
      {
        TrkrDefs::hitsetkey hitset_key = iter->first;
        TrkrHitSet* hitset = iter->second;

        // get all hits for this hitset
        TrkrHitSet::ConstRange hitrangei = hitset->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {
          TrkrDefs::hitkey hit_key = hitr->first;
          TrkrHit* hit = hitr->second;

          fx_hit[n_hit::nhitID] = hit_key;
          fx_hit[n_hit::nhite] = hit->getEnergy();
          fx_hit[n_hit::nhitadc] = hit->getAdc();
          unsigned int layer_local = TrkrDefs::getLayer(hitset_key);
          fx_hit[n_hit::nhitlayer] = (float) layer_local;
          fx_hit[n_hit::nhitphielem] = -666;
          fx_hit[n_hit::nhitzelem] = -666;

          if (layer_local >= 3 && layer_local < 7)
          {
            fx_hit[n_hit::nhitphielem] = InttDefs::getLadderPhiId(hitset_key);
            fx_hit[n_hit::nhitzelem] = InttDefs::getLadderZId(hitset_key);
          }
          if (layer_local >= 7 && layer_local < 55)
          {
            fx_hit[n_hit::nhitphielem] = TpcDefs::getSectorId(hitset_key);
            fx_hit[n_hit::nhitzelem] = TpcDefs::getSide(hitset_key);
          }
          /*
          if(layer_local>=55){
            if(MicromegasDefs::getSegmentationType(hitset_key)==MicromegasDefs::SEGMENTATION_Z){
              sector = 1;
              side = MicromegasDefs::getStrip(hit_key);
            }else{
              sector =MicromegasDefs::getStrip(hit_key);
              side =  1;
            }
          }
          */
          fx_hit[n_hit::nhitphielem] = TpcDefs::getSectorId(hitset_key);
          fx_hit[n_hit::nhitzelem] = TpcDefs::getSide(hitset_key);
          fx_hit[n_hit::nhitcellID] = 0;
          fx_hit[n_hit::nhitecell] = hit->getAdc();
          fx_hit[n_hit::nhitphibin] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhittbin] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhitphi] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhitr] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhitx] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhity] = std::numeric_limits<float>::quiet_NaN();
          fx_hit[n_hit::nhitz] = std::numeric_limits<float>::quiet_NaN();

          if (layer_local >= _nlayers_maps + _nlayers_intt && layer_local < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
          {
            PHG4TpcCylinderGeom* GeoLayer_local = _geom_container->GetLayerCellGeom(layer_local);
            double radius = GeoLayer_local->get_radius();
            fx_hit[n_hit::nhitphibin] = (float) TpcDefs::getPad(hit_key);
            fx_hit[n_hit::nhittbin] = (float) TpcDefs::getTBin(hit_key);
            fx_hit[n_hit::nhitphi] = GeoLayer_local->get_phicenter(fx_hit[n_hit::nhitphibin]);
            float phi = GeoLayer_local->get_phicenter(TpcDefs::getPad(hit_key));
            float clockperiod = GeoLayer_local->get_zstep();
            auto glob = m_tGeometry->getGlobalPositionTpc(hitset_key, hit_key, phi, radius, clockperiod);

            fx_hit[n_hit::nhitz] = glob.z();
            fx_hit[n_hit::nhitr] = radius;
            fx_hit[n_hit::nhitx] = glob.x();
            fx_hit[n_hit::nhity] = glob.y();
          }

          std::copy(fx_event, fx_event + n_event::evsize, hit_data);
          std::copy(fx_hit, fx_hit + n_hit::hitsize, hit_data + n_event::evsize);
          std::copy(fx_info, fx_info + ((int) (n_info::infosize)), hit_data + n_event::evsize + n_hit::hitsize);
          _ntp_hit->Fill(hit_data);
          // erase hit data
          std::copy(fx_hit_0, fx_hit_0 + ((int) (n_info::infosize)) + n_event::evsize + n_hit::hitsize, hit_data);
        }
      }
      delete[] hit_data;
    }
    if (Verbosity() >= 1)
    {
      _timer->stop();
      std::cout << "hit time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (Verbosity() >= 1)
  {
    std::cout << "check for ntp_cluster" << std::endl;
    _timer->restart();
  }

  if (_ntp_cluster)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_cluster (all of them) " << std::endl;
    }

    TrkrHitSetContainer* hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    //    TrkrClusterIterationMapv1* _iteration_map = findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    ClusterErrorPara ClusErrPara;

    if (Verbosity() > 1)
    {
      if (_cluster_map != nullptr)
      {
        std::cout << "got clustermap" << std::endl;
      }
      else
      {
        std::cout << "no clustermap" << std::endl;
      }

      if (hitsets != nullptr)
      {
        std::cout << "got hitsets" << std::endl;
      }
      else
      {
        std::cout << "no hitsets" << std::endl;
      }
    }

    if (_cluster_map && hitsets)
    {
      float* cluster_data = new float[((int) (n_info::infosize)) + n_event::evsize + n_cluster::clusize];
      Float_t fx_cluster_0[((int) (n_info::infosize)) + n_event::evsize + n_cluster::clusize] = {0};
      for (const auto& hitsetkey : _cluster_map->getHitSetKeys())
      {
        auto range = _cluster_map->getClusters(hitsetkey);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          TrkrDefs::cluskey cluster_key = iter->first;
          Float_t fx_cluster[n_cluster::clusize];
          FillCluster(&fx_cluster[0], cluster_key);

          std::copy(fx_event, fx_event + n_event::evsize, cluster_data);
          std::copy(fx_cluster, fx_cluster + n_cluster::clusize, cluster_data + n_event::evsize);
          std::copy(fx_info, fx_info + ((int) (n_info::infosize)), cluster_data + n_event::evsize + n_cluster::clusize);
          _ntp_cluster->Fill(cluster_data);
          std::copy(fx_cluster_0, fx_cluster_0 + ((int) (n_info::infosize)) + n_event::evsize + n_cluster::clusize, cluster_data);
        }
      }
      delete[] cluster_data;
    }
  }

  if (Verbosity() >= 1)
  {
    _timer->stop();
    std::cout << "cluster time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
  }
  if (_ntp_clus_trk)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_clus_trk " << std::endl;
      if (_cluster_map != nullptr)
      {
        std::cout << "got clustermap" << std::endl;
      }
      else
      {
        std::cout << "no clustermap" << std::endl;
      }
    }

    TrackSeedContainer* _tpc_seeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
    if (!_tpc_seeds)
    {
      std::cout << PHWHERE << " ERROR: Can't find "
                << "TpcTrackSeedContainer" << std::endl;
      return;
    }

    if (_trackmap)
    {
      int trackID = 0;
      float* clus_trk_data = new float[((int) (n_info::infosize)) + n_cluster::clusize + n_residual::ressize + n_seed::seedsize + n_event::evsize];
      float fx_clus_trk_0[((int) (n_info::infosize)) + n_cluster::clusize + n_residual::ressize + n_seed::seedsize + n_event::evsize] = {0};
      for (auto& iter : *_trackmap)
      {
        trackID++;
        SvtxTrack* track = iter.second;
        TrackSeed* tpcseed = track->get_tpc_seed();
        TrackSeed* siseed = track->get_silicon_seed();
        std::vector<Acts::Vector3> clusterPositions;
        std::vector<TrkrDefs::cluskey> clusterKeys;

        TrackFitUtils::position_vector_t xypoints;
        TrackFitUtils::position_vector_t rzpoints;

        clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                           tpcseed->end_cluster_keys());

        TrackFitUtils::getTrackletClusters(_tgeometry, _cluster_map,
                                           clusterPositions, clusterKeys);
        for (auto& pos : clusterPositions)
        {
          float clusr = sqrt((pos.x() * pos.x()) + (pos.y() * pos.y()));
          if (pos.y() < 0)
          {
            clusr *= -1;
          }

          // exclude silicon and tpot clusters for now
          if (std::fabs(clusr) > 80 || std::fabs(clusr) < 30)
          {
            continue;
          }
          rzpoints.emplace_back(pos.z(), clusr);
          xypoints.emplace_back(pos.x(), pos.y());
        }

        std::vector<float> fitparams = TrackFitUtils::fitClusters(clusterPositions, clusterKeys);
        if (fitparams.empty())
        {
          std::cout << "fit failed bailing...." << std::endl;
          continue;
        }

        float charge = std::numeric_limits<float>::quiet_NaN();
        if (tpcseed->get_qOverR() > 0)
        {
          charge = 1;
        }
        else
        {
          charge = -1;
        }

        //	      "pt:eta:phi:X0:Y0:charge:nhits:"
        float tpt = tpcseed->get_pt();
        float tptot = tpcseed->get_p();
        float teta = tpcseed->get_eta();
        float tphi = tpcseed->get_phi();
        auto xyparams = TrackFitUtils::line_fit(xypoints);
        auto rzparams = TrackFitUtils::line_fit(rzpoints);
        float xyint = std::get<1>(xyparams);
        float xyslope = std::get<0>(xyparams);
        float rzint = std::get<1>(rzparams);
        float rzslope = std::get<0>(rzparams);
        float R0 = abs(-1 * xyint) / std::sqrt((xyslope * xyslope) + 1);
        float tX0 = tpcseed->get_X0();
        float tY0 = tpcseed->get_Y0();
        float tZ0 = tpcseed->get_Z0();

        float nsil_local = 0;
        float ntpc_local = 0;
        if (siseed)
        {
          nsil_local = siseed->size_cluster_keys();
        }
        ntpc_local = tpcseed->size_cluster_keys();
        float nhits_local = clusterPositions.size();
        if (Verbosity() > 1)
        {
          std::cout << " tpc: " << tpcseed->size_cluster_keys() << std::endl;
          if (siseed)
          {
            std::cout << " si " << siseed->size_cluster_keys() << std::endl;
          }
        }
        //      nhits_local += tpcseed->size_cluster_keys();
        // fill the Gseed NTuple
        //---------------------
        float layerThicknesses[4] = {0.0, 0.0, 0.0, 0.0};
        // These are randomly chosen layer thicknesses for the TPC, to get the
        // correct region thicknesses in an easy to pass way to the helper fxn
        layerThicknesses[0] = _geom_container->GetLayerCellGeom(7)->get_thickness();
        layerThicknesses[1] = _geom_container->GetLayerCellGeom(8)->get_thickness();
        layerThicknesses[2] = _geom_container->GetLayerCellGeom(27)->get_thickness();
        layerThicknesses[3] = _geom_container->GetLayerCellGeom(50)->get_thickness();

        float dedx = TrackAnalysisUtils::calc_dedx(tpcseed, _cluster_map, _tgeometry, layerThicknesses);

        float pidedx = 0;
        float kdedx = 0;
        float prdedx = 0;
        if (charge > 0)
        {
          pidedx = f_pion_plus->Eval(tptot);
          kdedx = f_kaon_plus->Eval(tptot);
          prdedx = f_proton_plus->Eval(tptot);
        }
        else
        {
          pidedx = f_pion_minus->Eval(tptot);
          kdedx = f_kaon_minus->Eval(tptot);
          prdedx = f_proton_plus->Eval(tptot);
        }
        float n1pix = get_n1pix(tpcseed);
        float fx_seed[n_seed::seedsize] = {(float) trackID, 0, tpt, tptot, teta, tphi, xyint, rzint, xyslope, rzslope, tX0, tY0, tZ0, R0, charge, dedx, pidedx, kdedx, prdedx, n1pix, nsil_local, ntpc_local, nhits_local};
        if (_ntp_tpcseed)
        {
          float* tpcseed_data = new float[((int) (n_info::infosize)) + n_cluster::clusize + n_residual::ressize + n_seed::seedsize + n_event::evsize];
          std::copy(fx_event, fx_event + n_event::evsize, tpcseed_data);
          std::copy(fx_seed, fx_seed + n_seed::seedsize, tpcseed_data + n_event::evsize);
          std::copy(fx_info, fx_info + ((int) (n_info::infosize)), tpcseed_data + n_event::evsize + n_seed::seedsize);
          _ntp_tpcseed->Fill(tpcseed_data);
          delete[] tpcseed_data;
        }

        for (unsigned int i = 0; i < clusterPositions.size(); i++)
        {
          TrkrDefs::cluskey cluster_key = clusterKeys.at(i);
          unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
          Acts::Vector3 position = clusterPositions[i];
          Acts::Vector3 pca = TrackFitUtils::get_helix_pca(fitparams, position);
          float cluster_phi = atan2(position(1), position(0));
          float pca_phi = atan2(pca(1), pca(0));
          float dphi = cluster_phi - pca_phi;
          if (dphi > M_PI)
          {
            dphi = 2 * M_PI - dphi;
          }
          if (dphi < -M_PI)
          {
            dphi = 2 * M_PI + dphi;
          }

          float dz = position(2) - pca(2);

          float resr = sqrt((position(0) * position(0)) + (position(1) * position(1)) + (position(2) * position(2)));
          float seedR = std::abs(1.0 / tpcseed->get_qOverR());
          float alpha = (resr * resr) / (2 * resr * seedR);
          float beta = std::abs(std::atan(tpcseed->get_slope()));
          float fx_cluster[n_cluster::clusize];
          if (layer_local >= 7 && layer_local < 55)
          {
            PHG4TpcCylinderGeom* GeoLayer_local = _geom_container->GetLayerCellGeom(layer_local);
            float thick = GeoLayer_local->get_thickness();

            float alphacorr = std::cos(alpha);
            if (alphacorr < 0 || alphacorr > 4)
            {
              alphacorr = 4;
            }
            float betacorr = std::cos(beta);
            if (betacorr < 0 || betacorr > 4)
            {
              betacorr = 4;
            }
            fx_cluster[ncluthick] = thick;
            fx_cluster[ncluafac] = alphacorr;
            fx_cluster[nclubfac] = betacorr;
          }
          else
          {
            fx_cluster[ncluthick] = 0;
            fx_cluster[ncluafac] = 0;
            fx_cluster[nclubfac] = 0;
          }
          float fx_res[n_residual::ressize] = {alpha, beta, dphi, dphi, dz};
          // sphi:syxint:srzint:sxyslope:srzslope:sX0:sY0:sdZ0:sR0
          FillCluster(&fx_cluster[0], cluster_key);
          std::copy(fx_event, fx_event + n_event::evsize, clus_trk_data);
          std::copy(fx_cluster, fx_cluster + n_cluster::clusize, clus_trk_data + n_event::evsize);
          std::copy(fx_res, fx_res + n_residual::ressize, clus_trk_data + n_event::evsize + n_cluster::clusize);
          std::copy(fx_seed, fx_seed + n_seed::seedsize, clus_trk_data + n_event::evsize + n_cluster::clusize + n_residual::ressize);
          std::copy(fx_info, fx_info + ((int) (n_info::infosize)), clus_trk_data + n_event::evsize + n_cluster::clusize + n_residual::ressize + n_seed::seedsize);
          _ntp_clus_trk->Fill(clus_trk_data);
          std::copy(fx_clus_trk_0, fx_clus_trk_0 + ((int) (n_info::infosize)) + n_cluster::clusize + n_residual::ressize + n_seed::seedsize + n_event::evsize, clus_trk_data);
        }
      }
      delete[] clus_trk_data;
    }
  }

  //------------------------
  // fill the Track NTuple
  //------------------------

  if (_ntp_track)
  {
    if (Verbosity() >= 1)
    {
      std::cout << "Filling ntp_track " << std::endl;
      _timer->restart();
    }
    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (_trackmap)
    {
      for (auto& iter : *_trackmap)
      {
        SvtxTrack* track = iter.second;
        float fx_track[n_track::trksize];
        FillTrack(&fx_track[0], track, vertexmap);
        float* track_data = new float[((int) (n_info::infosize)) + n_track::trksize + n_event::evsize];
        std::copy(fx_event, fx_event + n_event::evsize, track_data);
        std::copy(fx_track, fx_track + n_track::trksize, track_data + n_event::evsize);
        std::copy(fx_info, fx_info + ((int) (n_info::infosize)), track_data + n_event::evsize + n_track::trksize);
        _ntp_track->Fill(track_data);
        delete[] track_data;
      }
    }
    if (Verbosity() > 1)
    {
      _timer->stop();
      std::cout << "track time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }

  if (_ntp_siseed)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Filling ntp_tpcseed " << std::endl;
      _timer->restart();
    }

    if (Verbosity() > 1)
    {
      _timer->stop();
      std::cout << "tpcseed time:                " << _timer->get_accumulated_time() / 1000. << " sec" << std::endl;
    }
  }
  return;
}

void TrkrNtuplizer::FillTrack(float fX[50], SvtxTrack* track, GlobalVertexMap* vertexmap)
{
  float trackID = track->get_id();
  TrackSeed* tpcseed = track->get_tpc_seed();
  TrackSeed* silseed = track->get_silicon_seed();
  short int crossing_int = track->get_crossing();
  fX[n_track::ntrktrackID] = trackID;
  fX[n_track::ntrkcrossing] = std::numeric_limits<float>::quiet_NaN();
  if (crossing_int == SHRT_MAX)
  {
    fX[n_track::ntrkcrossing] = std::numeric_limits<float>::quiet_NaN();
  }
  else
  {
    fX[n_track::ntrkcrossing] = (float) crossing_int;
  }
  fX[n_track::ntrkcharge] = track->get_charge();
  fX[n_track::ntrkquality] = track->get_quality();
  fX[n_track::ntrkchisq] = track->get_chisq();
  fX[n_track::ntrkndf] = track->get_ndf();
  fX[n_track::ntrknhits] = 0;
  if (tpcseed)
  {
    fX[n_track::ntrknhits] += tpcseed->size_cluster_keys();
  }
  if (silseed)
  {
    fX[n_track::ntrknhits] += silseed->size_cluster_keys();
  }

  fX[n_track::ntrknmaps] = 0;
  fX[n_track::ntrknintt] = 0;
  fX[n_track::ntrknmms] = 0;
  fX[n_track::ntrkntpc] = 0;
  fX[n_track::ntrkntpc1] = 0;
  fX[n_track::ntrkntpc11] = 0;
  fX[n_track::ntrkntpc2] = 0;
  fX[n_track::ntrkntpc3] = 0;
  fX[n_track::ntrkndedx] = 0;
  fX[n_track::ntrknpidedx] = 0;
  fX[n_track::ntrknkdedx] = 0;
  fX[n_track::ntrknprdedx] = 0;
  if (tpcseed)
  {
    float layerThicknesses[4] = {0.0, 0.0, 0.0, 0.0};
    // These are randomly chosen layer thicknesses for the TPC, to get the
    // correct region thicknesses in an easy to pass way to the helper fxn
    layerThicknesses[0] = _geom_container->GetLayerCellGeom(7)->get_thickness();
    layerThicknesses[1] = _geom_container->GetLayerCellGeom(8)->get_thickness();
    layerThicknesses[2] = _geom_container->GetLayerCellGeom(27)->get_thickness();
    layerThicknesses[3] = _geom_container->GetLayerCellGeom(50)->get_thickness();

    fX[n_track::ntrkndedx] = TrackAnalysisUtils::calc_dedx(tpcseed, _cluster_map, _tgeometry, layerThicknesses);
    
    float trptot = track->get_p();
    if (track->get_charge() > 0)
    {
      fX[n_track::ntrknpidedx] = f_pion_plus->Eval(trptot);
      fX[n_track::ntrknkdedx] = f_kaon_plus->Eval(trptot);
      fX[n_track::ntrknprdedx] = f_proton_plus->Eval(trptot);
    }
    else
    {
      fX[n_track::ntrknpidedx] = f_pion_minus->Eval(trptot);
      fX[n_track::ntrknkdedx] = f_kaon_minus->Eval(trptot);
      fX[n_track::ntrknprdedx] = f_proton_minus->Eval(trptot);
    }

    for (SvtxTrack::ConstClusterKeyIter iter_local = tpcseed->begin_cluster_keys();
         iter_local != tpcseed->end_cluster_keys();
         ++iter_local)
    {
      TrkrDefs::cluskey cluster_key = *iter_local;
      unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
      switch (TrkrDefs::getTrkrId(cluster_key))
      {
      case TrkrDefs::TrkrId::mvtxId:
        fX[n_track::ntrknmaps]++;
        break;
      case TrkrDefs::TrkrId::inttId:
        fX[n_track::ntrknintt]++;
        break;
      case TrkrDefs::TrkrId::tpcId:
        fX[n_track::ntrkntpc]++;
        if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 8)
        {
          fX[n_track::ntrkntpc11]++;
        }
        if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 16)
        {
          fX[n_track::ntrkntpc1]++;
        }
        else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 32)
        {
          fX[n_track::ntrkntpc2]++;
        }
        else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 48)
        {
          fX[n_track::ntrkntpc3]++;
        }
        break;
      case TrkrDefs::TrkrId::micromegasId:
        fX[n_track::ntrknmms]++;
        break;
      default:
        break;
      }
    }
  }
  if (silseed)
  {
    for (SvtxTrack::ConstClusterKeyIter iter_local = silseed->begin_cluster_keys();
         iter_local != silseed->end_cluster_keys();
         ++iter_local)
    {
      TrkrDefs::cluskey cluster_key = *iter_local;
      unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
      switch (TrkrDefs::getTrkrId(cluster_key))
      {
      case TrkrDefs::TrkrId::mvtxId:
        fX[n_track::ntrknmaps]++;
        break;
      case TrkrDefs::TrkrId::inttId:
        fX[n_track::ntrknintt]++;
        break;
      case TrkrDefs::TrkrId::tpcId:
        fX[n_track::ntrkntpc]++;
        if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 8)
        {
          fX[n_track::ntrkntpc11]++;
        }
        if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 16)
        {
          fX[n_track::ntrkntpc1]++;
        }
        else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 32)
        {
          fX[n_track::ntrkntpc2]++;
        }
        else if ((layer_local - (_nlayers_maps + _nlayers_intt)) < 48)
        {
          fX[n_track::ntrkntpc3]++;
        }
        break;
      case TrkrDefs::TrkrId::micromegasId:
        fX[n_track::ntrknmms]++;
        break;
      default:
        break;
      }
    }
  }
  fX[n_track::ntrkdca3dxy] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca3dz] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca3dxysigma] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca3dzsigma] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca2d] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkdca2dsigma] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkvx] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkvy] = std::numeric_limits<float>::quiet_NaN();
  fX[n_track::ntrkvz] = std::numeric_limits<float>::quiet_NaN();

  int vertexID = track->get_vertex_id();
  fX[n_track::ntrkvertexID] = vertexID;

  if (vertexID >= 0 && vertexmap != nullptr)
  {
    if (vertexmap->size() > 100000)
    {
      std::cout << "too many vtx's" << std::endl;
    }
    /*
    auto vertexit = vertexmap->find(vertexID);
    if (vertexit != vertexmap->end())
      {
        auto vertex = vertexit->second;

        float vx = vertex->get_x();
        float vy = vertex->get_y();
        float vz = vertex->get_z();
        fX[n_track::ntrkvx] = vx;
        fX[n_track::ntrkvy] = vy;
        fX[n_track::ntrkvz] = vz;
        Acts::Vector3 vert(vx, vy, vz);
        Acts::Vector3 zero = Acts::Vector3::Zero();
        auto dcapair = TrackAnalysisUtils::get_dca(track, zero);
        fX[n_track::ntrkdca3dxy] = dcapair.first.first;
        fX[n_track::ntrkdca3dxysigma] = dcapair.first.second;
        fX[n_track::ntrkdca3dz] = dcapair.second.first;
        fX[n_track::ntrkdca3dzsigma] = dcapair.second.second;

      }
    */
  }

  fX[n_track::ntrkcharge] = track->get_charge();
  float px = track->get_px();
  float py = track->get_py();
  float pz = track->get_pz();
  fX[n_track::ntrkpx] = px;
  fX[n_track::ntrkpy] = py;
  fX[n_track::ntrkpz] = pz;
  TVector3 v(px, py, pz);
  fX[n_track::ntrkpt] = v.Pt();
  fX[n_track::ntrketa] = v.Eta();
  fX[n_track::ntrkphi] = v.Phi();
  float CVxx = track->get_error(3, 3);
  float CVxy = track->get_error(3, 4);
  float CVxz = track->get_error(3, 5);
  float CVyy = track->get_error(4, 4);
  float CVyz = track->get_error(4, 5);
  float CVzz = track->get_error(5, 5);
  fX[n_track::ntrkdeltapt] = std::sqrt((CVxx * px * px + 2 * CVxy * px * py + CVyy * py * py) / (px * px + py * py));
  fX[n_track::ntrkdeltaeta] = std::sqrt((CVzz * (px * px + py * py) * (px * px + py * py) + pz * (-2 * (CVxz * px + CVyz * py) * (px * px + py * py) + CVxx * px * px * pz + CVyy * py * py * pz + 2 * CVxy * px * py * pz)) / ((px * px + py * py) * (px * px + py * py) * (px * px + py * py + pz * pz)));
  fX[n_track::ntrkdeltaphi] = std::sqrt((CVyy * px * px - 2 * CVxy * px * py + CVxx * py * py) / ((px * px + py * py) * (px * px + py * py)));

  fX[n_track::ntrkpcax] = track->get_x();
  fX[n_track::ntrkpcay] = track->get_y();
  fX[n_track::ntrkpcaz] = track->get_z();
}

float TrkrNtuplizer::get_n1pix(TrackSeed* tpcseed)
{
  std::vector<TrkrDefs::cluskey> clusterKeys;
  clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                     tpcseed->end_cluster_keys());

  float n1pix = 0;
  for (unsigned long cluster_key : clusterKeys)
  {
    TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);
    if (cluster->getPhiSize() == 1)
    {
      n1pix++;
    }
  }
  return n1pix;
}

void TrkrNtuplizer::FillCluster(float fXcluster[n_cluster::clusize], TrkrDefs::cluskey cluster_key)
{
  unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
  TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);

  Acts::Vector3 cglob;
  cglob = _tgeometry->getGlobalPosition(cluster_key, cluster);
  float x = cglob(0);
  float y = cglob(1);
  float z = cglob(2);
  TVector3 pos(x, y, z);
  float r = pos.Perp();
  float phi = pos.Phi();
  auto para_errors = ClusterErrorPara::get_clusterv5_modified_error(cluster, r, cluster_key);
  unsigned int phibin = std::numeric_limits<unsigned int>::quiet_NaN();
  float tbin = std::numeric_limits<float>::quiet_NaN();
  float locx = cluster->getLocalX();
  float locy = cluster->getLocalY();
  fXcluster[n_cluster::ncludcal] = 1;
  if (layer_local >= _nlayers_maps + _nlayers_intt && layer_local < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
  {
    int side_tpc = TpcDefs::getSide(cluster_key);
    PHG4TpcCylinderGeom* GeoLayer_local = _geom_container->GetLayerCellGeom(layer_local);
    // NOLINTNEXTLINE(bugprone-incorrect-roundings)
    phibin = (unsigned int) GeoLayer_local->get_pad_float(phi, side_tpc) + 0.5;
    tbin = GeoLayer_local->get_tbin_float(locy - 39.6);
    unsigned int cside = TpcDefs::getSide(cluster_key);
    unsigned int csector = TpcDefs::getSectorId(cluster_key);
    unsigned int csegment = 0;
    TrkrDefs::cluskey dummykey = TpcDefs::genClusKey(layer_local, csector, cside, phibin);
    std::map<TrkrDefs::cluskey, fee_info>::iterator fee_it = fee_map.find(dummykey);
    if (fee_it != fee_map.end())
    {
      fXcluster[n_cluster::nclufee] = (*fee_it).second.fee;
      fXcluster[n_cluster::ncluchan] = (*fee_it).second.channel;
      fXcluster[n_cluster::nclusampa] = (*fee_it).second.sampa;
    }
    else
    {
      fXcluster[n_cluster::nclufee] = 0;
      fXcluster[n_cluster::ncluchan] = 0;
    }
    if (layer_local >= 7 + 32)
    {
      csegment = 2;
    }
    else if (layer_local >= 7 + 16)
    {
      csegment = 1;
    }
    if (_do_dedx_calib == true)
    {
      fXcluster[n_cluster::ncludcal] = dedxcorr[cside][csector][csegment];
    }
    else
    {
      fXcluster[n_cluster::ncludcal] = 1;
    }
  }
  else
  {
    phibin = locx;
    tbin = locy;
  }
  fXcluster[n_cluster::nclulocx] = locx;
  fXcluster[n_cluster::nclulocy] = locy;
  fXcluster[n_cluster::nclux] = cglob(0);
  fXcluster[n_cluster::ncluy] = cglob(1);
  fXcluster[n_cluster::ncluz] = cglob(2);
  fXcluster[n_cluster::nclur] = pos.Perp();
  fXcluster[n_cluster::ncluphi] = pos.Phi();
  fXcluster[n_cluster::nclueta] = pos.Eta();
  fXcluster[n_cluster::nclutheta] = pos.Theta();

  fXcluster[n_cluster::ncluphibin] = phibin;
  fXcluster[n_cluster::nclutbin] = tbin;

  fXcluster[n_cluster::ncluex] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::ncluey] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::ncluez] = sqrt(para_errors.second);
  fXcluster[n_cluster::ncluephi] = sqrt(para_errors.first);
  fXcluster[n_cluster::nclupez] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::nclupephi] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::nclue] = cluster->getAdc();
  fXcluster[n_cluster::ncluadc] = cluster->getAdc();
  fXcluster[n_cluster::nclumaxadc] = cluster->getMaxAdc();
  fXcluster[n_cluster::nclulayer] = layer_local;

  if (layer_local < 3)
  {
    fXcluster[n_cluster::ncluphielem] = MvtxDefs::getStaveId(cluster_key);
    fXcluster[n_cluster::ncluzelem] = MvtxDefs::getChipId(cluster_key);
  }
  if (layer_local >= 3 && layer_local < 7)
  {
    fXcluster[n_cluster::ncluphielem] = InttDefs::getLadderPhiId(cluster_key);
    fXcluster[n_cluster::ncluzelem] = InttDefs::getLadderZId(cluster_key);
  }
  if (layer_local >= 7 && layer_local < 55)
  {
    fXcluster[n_cluster::ncluphielem] = TpcDefs::getSectorId(cluster_key);
    fXcluster[n_cluster::ncluzelem] = TpcDefs::getSide(cluster_key);
  } /*
      if(layer_local>=55){
      if(MicromegasDefs::getSegmentationType(hitsetkey)==MicromegasDefs::SEGMENTATION_Z){
      sector = 1;
      side = MicromegasDefs::getTileId(hitsetkey);
      }else{
      sector =MicromegasDefs::getTileId(hitsetkey);
      side =  1;
      }
      }
    */
  fXcluster[n_cluster::nclusize] = cluster->getSize();
  fXcluster[n_cluster::ncluphisize] = cluster->getPhiSize();
  fXcluster[n_cluster::ncluzsize] = cluster->getZSize();
  fXcluster[n_cluster::nclupedge] = cluster->getEdge();
  if (layer_local == 7 || layer_local == 22 ||
      layer_local == 23 || layer_local == 28 ||
      layer_local == 39 || layer_local == 54)
  {
    fXcluster[n_cluster::ncluredge] = 1;
  }

  fXcluster[n_cluster::ncluovlp] = 3;  // cluster->getOvlp();
  fXcluster[n_cluster::nclutrackID] = std::numeric_limits<float>::quiet_NaN();
  fXcluster[n_cluster::ncluniter] = 0;

  return;
}

std::vector<TrkrDefs::cluskey> TrkrNtuplizer::get_track_ckeys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> cluster_keys;
  TrackSeed* tpcseed = track->get_tpc_seed();
  TrackSeed* silseed = track->get_silicon_seed();
  if (silseed)
  {
    for (auto iter = silseed->begin_cluster_keys();
         iter != silseed->end_cluster_keys();
         ++iter)
    {
      cluster_keys.push_back(*iter);
    }
  }
  if (tpcseed)
  {
    for (auto iter = tpcseed->begin_cluster_keys();
         iter != tpcseed->end_cluster_keys();
         ++iter)
    {
      cluster_keys.push_back(*iter);
    }
  }

  return cluster_keys;
}

SvtxTrack* TrkrNtuplizer::best_track_from(TrkrDefs::cluskey cluster_key)
{
  std::map<TrkrDefs::cluskey, SvtxTrack*>::iterator find_iter =
      _cache_best_track_from_cluster.find(cluster_key);
  if (find_iter != _cache_best_track_from_cluster.end())
  {
    return find_iter->second;
  }

  SvtxTrack* best_track = nullptr;
  float best_quality = FLT_MAX;

  std::set<SvtxTrack*> tracks = all_tracks_from(cluster_key);
  // loop over all SvtxTracks
  for (auto* candidate : tracks)
  {
    if (candidate->get_quality() < best_quality)
    {
      best_quality = candidate->get_quality();
      best_track = candidate;
    }
  }

  _cache_best_track_from_cluster.insert(std::make_pair(cluster_key, best_track));
  return best_track;
}

std::set<SvtxTrack*> TrkrNtuplizer::all_tracks_from(TrkrDefs::cluskey cluster_key)
{
  std::set<SvtxTrack*> tracks;

  if (_cache_track_from_cluster_exists == false)
  {
    create_cache_track_from_cluster();
  }
  std::map<TrkrDefs::cluskey, std::set<SvtxTrack*> >::iterator find_iter =
      _cache_all_tracks_from_cluster.find(cluster_key);
  if (find_iter != _cache_all_tracks_from_cluster.end())
  {
    return find_iter->second;
  }

  return tracks;

  // loop over all SvtxTracks
  for (auto& iter : *_trackmap)
  {
    SvtxTrack* track = iter.second;
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);

    // loop over all clusters
    for (const auto& candidate : cluster_keys)
    {
      //      if (_strict)
      //      {
      //        assert(candidate);
      //      }
      //      else if (!candidate)
      //      {
      //        ++_errors;
      //        continue;
      //      }

      if (cluster_key == candidate)
      {
        tracks.insert(track);
      }
    }
  }

  _cache_all_tracks_from_cluster.insert(std::make_pair(cluster_key, tracks));

  return tracks;
}

void TrkrNtuplizer::create_cache_track_from_cluster()
{
  if (!_trackmap)
  {
    return;
  }

  // loop over all SvtxTracks
  for (auto& iter : *_trackmap)
  {
    SvtxTrack* track = iter.second;
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);

    // loop over all clusters
    for (const auto& candidate_key : cluster_keys)
    {
      // check if cluster has an entry in cache
      std::map<TrkrDefs::cluskey, std::set<SvtxTrack*> >::iterator cliter =
          _cache_all_tracks_from_cluster.find(candidate_key);
      if (cliter != _cache_all_tracks_from_cluster.end())
      {                                // got entry
        cliter->second.insert(track);  // add track to list;
      }
      else
      {
        std::set<SvtxTrack*> tracks;
        tracks.insert(track);
        _cache_all_tracks_from_cluster.insert(std::make_pair(candidate_key, tracks));
      }
    }
  }
  _cache_track_from_cluster_exists = true;

  return;
}

TMatrixF TrkrNtuplizer::calculateClusterError(TrkrCluster* c, float& clusphi)
{
  TMatrixF localErr(3, 3);
  localErr[0][0] = 0.;
  localErr[0][1] = 0.;
  localErr[0][2] = 0.;
  localErr[1][0] = 0.;
  localErr[1][1] = c->getActsLocalError(0, 0);
  localErr[1][2] = c->getActsLocalError(0, 1);
  localErr[2][0] = 0.;
  localErr[2][1] = c->getActsLocalError(1, 0);
  localErr[2][2] = c->getActsLocalError(1, 1);

  TMatrixF ROT(3, 3);
  ROT[0][0] = std::cos(clusphi);
  ROT[0][1] = -std::sin(clusphi);
  ROT[0][2] = 0.0;
  ROT[1][0] = std::sin(clusphi);
  ROT[1][1] = std::cos(clusphi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0;
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;
  TMatrixF ROT_T(3, 3);
  ROT_T.Transpose(ROT);

  TMatrixF err(3, 3);
  err = ROT * localErr * ROT_T;
  return err;
}
