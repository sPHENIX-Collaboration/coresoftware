#include "JetEvaluator.h"

#include "JetEvalStack.h"
#include "JetRecoEval.h"

#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TNtuple.h>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>

using namespace std;

JetEvaluator::JetEvaluator(const string &name,
                           const string &recojetname,
                           const string &truthjetname,
                           const string &filename)
  : SubsysReco(name)
  , _recojetname(recojetname)
  , _truthjetname(truthjetname)
  , _ievent(0)
  , _jetevalstack(nullptr)
  , _strict(false)
  , _do_recojet_eval(true)
  , _do_truthjet_eval(true)
  , _ntp_recojet(nullptr)
  , _ntp_truthjet(nullptr)
  , _filename(filename)
  , _tfile(nullptr)
{
}

int JetEvaluator::Init(PHCompositeNode */*topNode*/)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

  if (_do_recojet_eval) _ntp_recojet = new TNtuple("ntp_recojet", "reco jet => max truth jet",
                                                   "event:id:ncomp:eta:phi:e:pt:"
                                                   "gid:gncomp:geta:gphi:ge:gpt:"
                                                   "efromtruth");

  if (_do_truthjet_eval) _ntp_truthjet = new TNtuple("ntp_truthjet", "truth jet => best reco jet",
                                                     "event:gid:gncomp:geta:gphi:ge:gpt:"
                                                     "id:ncomp:eta:phi:e:pt:"
                                                     "efromtruth");

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetEvaluator::process_event(PHCompositeNode *topNode)
{
  if (!_jetevalstack)
  {
    _jetevalstack = new JetEvalStack(topNode, _recojetname, _truthjetname);
    _jetevalstack->set_strict(_strict);
    _jetevalstack->set_verbosity(Verbosity() + 1);
  }
  else
  {
    _jetevalstack->next_event(topNode);
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

  printOutputInfo(topNode);

  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetEvaluator::End(PHCompositeNode */*topNode*/)
{
  _tfile->cd();

  if (_do_recojet_eval) _ntp_recojet->Write();
  if (_do_truthjet_eval) _ntp_truthjet->Write();

  _tfile->Close();

  delete _tfile;

  if (Verbosity() > 0)
  {
    cout << "========================== JetEvaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  delete _jetevalstack;

  return Fun4AllReturnCodes::EVENT_OK;
}

void JetEvaluator::printInputInfo(PHCompositeNode */*topNode*/)
{
  // to be implemented later if needed
  return;
}

void JetEvaluator::printOutputInfo(PHCompositeNode */*topNode*/)
{
  // to be implemented later if needed
  return;
}

void JetEvaluator::fillOutputNtuples(PHCompositeNode *topNode)
{
  if (Verbosity() > 2) cout << "JetEvaluator::fillOutputNtuples() entered" << endl;

  JetRecoEval *recoeval = _jetevalstack->get_reco_eval();
  //JetTruthEval* trutheval = _jetevalstack->get_truth_eval();

  //-------------------------
  // fill the reco jet ntuple
  //-------------------------

  if (_do_recojet_eval)
  {
    if (Verbosity() > 1) cout << "JetEvaluator::filling recojet ntuple..." << endl;

    JetMap *recojets = findNode::getClass<JetMap>(topNode, _recojetname.c_str());
    if (!recojets)
    {
      cerr << PHWHERE << " ERROR: Can't find " << _recojetname << endl;
      exit(-1);
    }

    // for every recojet
    for (JetMap::Iter iter = recojets->begin();
         iter != recojets->end();
         ++iter)
    {
      Jet *recojet = iter->second;
      Jet *truthjet = recoeval->max_truth_jet_by_energy(recojet);

      float id = recojet->get_id();
      float ncomp = recojet->size_comp();
      float eta = recojet->get_eta();
      float phi = recojet->get_phi();
      float e = recojet->get_e();
      float pt = recojet->get_pt();

      float gid = NAN;
      float gncomp = NAN;
      float geta = NAN;
      float gphi = NAN;
      float ge = NAN;
      float gpt = NAN;
      float efromtruth = NAN;

      if (truthjet)
      {
        gid = truthjet->get_id();
        gncomp = truthjet->size_comp();
        geta = truthjet->get_eta();
        gphi = truthjet->get_phi();
        ge = truthjet->get_e();
        gpt = truthjet->get_pt();
        efromtruth = recoeval->get_energy_contribution(recojet, truthjet);
      }

      float recojet_data[14] = {(float) _ievent,
                                id,
                                ncomp,
                                eta,
                                phi,
                                e,
                                pt,
                                gid,
                                gncomp,
                                geta,
                                gphi,
                                ge,
                                gpt,
                                efromtruth};

      _ntp_recojet->Fill(recojet_data);
    }
  }

  //-------------------------
  // fill the truth jet ntuple
  //-------------------------

  if (_do_truthjet_eval)
  {
    if (Verbosity() > 1) cout << "JetEvaluator::filling truthjet ntuple..." << endl;

    JetMap *truthjets = findNode::getClass<JetMap>(topNode, _truthjetname.c_str());
    if (!truthjets)
    {
      cerr << PHWHERE << " ERROR: Can't find " << _truthjetname << endl;
      exit(-1);
    }

    // for every truthjet
    for (JetMap::Iter iter = truthjets->begin();
         iter != truthjets->end();
         ++iter)
    {
      Jet *truthjet = iter->second;
      Jet *recojet = recoeval->best_jet_from(truthjet);

      float gid = truthjet->get_id();
      float gncomp = truthjet->size_comp();
      float geta = truthjet->get_eta();
      float gphi = truthjet->get_phi();
      float ge = truthjet->get_e();
      float gpt = truthjet->get_pt();

      float id = NAN;
      float ncomp = NAN;
      float eta = NAN;
      float phi = NAN;
      float e = NAN;
      float pt = NAN;
      float efromtruth = NAN;

      if (recojet)
      {
        id = recojet->get_id();
        ncomp = recojet->size_comp();
        eta = recojet->get_eta();
        phi = recojet->get_phi();
        e = recojet->get_e();
        pt = recojet->get_pt();
        efromtruth = recoeval->get_energy_contribution(recojet, truthjet);
      }

      float truthjet_data[14] = {(float) _ievent,
                                 gid,
                                 gncomp,
                                 geta,
                                 gphi,
                                 ge,
                                 gpt,
                                 id,
                                 ncomp,
                                 eta,
                                 phi,
                                 e,
                                 pt,
                                 efromtruth};

      _ntp_truthjet->Fill(truthjet_data);
    }
  }

  return;
}
