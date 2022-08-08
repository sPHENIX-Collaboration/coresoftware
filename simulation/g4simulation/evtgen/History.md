# History file for EvtGen

From version 2.0.0,
Tabc labels refer to [Maniphest tasks](https://phab.hepforge.org/maniphest/query/nkBRd9OhPCBN/), while
Dxyz labels refer to [Differential code reviews](https://phab.hepforge.org/differential/query/YDY8EgjNGd.e/) on HepForge:

https://phab.hepforge.org/Tabc

https://phab.hepforge.org/Dxyz

===
## R02-01-01

8th Sep 2021 Michal Kreps
* Update location of web page for Pythia8 download in setup script.

8th Sep 2021 Markus Prim, Lu Cao, Chaoyi Lyu and Michel De Cian (Michal Kreps)
* D73: Add new model for semileptonic B decays with BCL and BGL form-factors

8th June 2021 Michal Kreps
* T110, D71: Fix B+ --> eta' l nu BF which was order of magnitude too high.  Balance the decrease by increasing B+ --> D0 l nu, which is after change still bit smaller than PDG 2021.

8th Jun 2021 Michal Kreps
* D71: Fix B+ --> eta' l nu BF.

20th Apr 2021 Tom Lathem
* D68: Fix compilation with Pythia 8.304

17th Mar 2021 Michal Kreps
* D62: Improve PI0DALITZ model to dynamically work out maximum probability to make it usuable also for eta --> llgamma decays. Model ETA2MUMUGAMMA is pure one-to-one copy of the PI0DALITZ and as such it is removed.

15th Jan 2021 Michal Kreps
* D47: Model to generate 4-body phase-space decays in restricted part of the m12-m34 space

12th Jan 2021 Michal Kreps
* D48: Fix bug in calculation of the daughter momentum in decay model EvtBsMuMuKK

7th Jan 2021 Michal Kreps
* D43: Allow to pass particle properties table in form of stringstream to constructor of EvtGen for use case where these are created on fly.

10th Dec 2020 Michal Kreps
* D36: EvtFlatSqDalitz model to be more efficient and to avoid cut-off around the edges

21st Aug 2020 John Back
* T109: Add EvtEtaLLPiPi model for eta' -> l l pi pi decays (l = e or mu),
  - courtesy of Aleksei Luchinsky (LHCb).

29th Jun 2020 Michal Kreps
* T103: Add missing include in EvtGenBase/EvtMatrix.hh.

15th May 2020 Michal Kreps
* D28: Add EvtThreeBodyPhsp (rectangular Dalitz plot selection) and
    EvtPhspDecaytimeCut (minimum decay time requirement) models.
* D27: Fix initialisation of constants in EvtBTo3hCP model.

===
## R02-00-00

24th Apr 2020 Michal Kreps
* Update particle properties file evt.pdl to 2019 version of RPP by PDG.

23rd Apr 2020 Tom Latham
* Apply copyright and licence notices to all relevant files.

17th Apr 2020 Tom Latham
* Add text of GNU GPLv3 in COPYING, add (preliminary) list of authors in
  AUTHORS, and clean-up all source files, prior to applying copyright and
  licence notices.

9th Apr 2020 Tom Latham
* Improve, and document use of, standalone installation script.
* Apply clang-format formatting to all C++ source files.

8th Apr 2020 Tom Latham
* One small modernisation change to EvtPhotosEngine to match that already
  applied in EvtTauolaEngine.

8th Apr 2020 John Back
* Fixed NonReson amplitude and the 4-momentum boosts used for angles in
  EvtLambdacPHH,
  - courtesy of Elisabeth Niel (LHCb).

7th Apr 2020 Gerhard Raven, Tom Latham, Michal Kreps and John Back
* Incorporate C++ modernization changes from Gerhard Raven (LHCb).
  - Merged modernize branch into master.

9th Mar 2020 John Back
* Add EvtAbsExternalGen::getDecayProb() to allow external generators to
  return a probability that can be used in EvtDecayProb (default = 1).

6th Mar 2020 Andrii Verbytskyi and Tom Latham
* Implement HepMC3 support: EvtHepMCEvent, external engines & cmake files.

15th Nov 2019 John Back
* Added EvtPsi2JpsiPiPi model for psi2S -> J/psi pi+ pi- decays based on hep-ph/1507.07985,
  - courtesy of Aleksei Luchinsky (LHCb).

21st August 2019 Michal Kreps
* Added the EvtDToKpienu decay model for D -> K pi e nu decays from BESIII,
  - courtesy of Liaoyuan Dong.

16th July 2019 John Back
* Correct imaginary sign for EvtComplex /= (EvtComplex c) operator.

3rd July 2019 John Back
* Added the EvtLambdacPHH decay model for Lc -> p K pi decays with K*(890), 
  Delta++(1232) and Lambda(1520) resonances, based on the Fermilab E791
  analysis hep-ex/9912003v1,
  - courtesy of Elisabeth Niel and Patrick Robbe (LHCb).
* Modified EvtResonance2 to accept other barrier factor radii.

3rd July 2019 Michal Kreps
* Make sure minimum mass for resonances with non-zero widths is larger than
  1e-4 GeV in EvtRelBreitWignerBarrierFact.

3rd May 2019 John Back
* Corrected EvtSLDiBaryonAmp bugs/issues in the BToDiBaryonlnupQCD model: 
  - parity, amplitude terms and B momentum reference frame variables.    
* Corrected treament of spinor indices in EvtRareLb2Lll,
  - courtesy of Tom Blake and Michal Kreps (LHCb).
* Updated the EvtBcVHad model to also handle Bc -> psi Ks K decays, 
  - courtesy of Aleksei Luchinsky (LHCb).
* Add new decay model EvtBsMuMuKK (BS_MUMUKK) for Bs to J/psi (mu+mu-) K+K-,
  - courtesy of Veronika Chobanova, Jeremy Dalseno, Diego Martinez Santos and
    Marcos Romero Lamas (LHCb).
* Fix infinite loop during initialisation of the EvtBTo3hCP model via
  EvtCBTo3piP00 and EvtCBTo3piMPP,
  - courtesy of Peter Richardson (Durham).

15th March 2019 Tom Latham
* Implement cmake build system, replacing the old config method.

30th Jan 2019 John Back
* Fix modernization compiler errors and warnings.

29th Jan 2019 Michal Kreps
* Allow reading decay files which are missing end-of-line before end-of-file.

21st December 2018 John Back
* Imported C++ modernization changes from Gerhard Raven (LHCb).

7th December 2018 John Back
* Added the EvtBLLNuL (BLLNUL) model that generates rare B -> ell ell nu ell 
  decays, where ell = e or mu,
  - courtesy of Anna Danilina and Nikolai Nikitin (LHCb).
* Removed the EvtB2MuMuMuNu (BUTOMMMN) model, since its now replaced 
  by the more general BLLNuL one.

5th November 2018 John Back
* Added the BToDiBaryonlnupQCD model for generating B to p N* l nu decays,
  where N can be any (exited) charged baryon (spin 1/2 or 3/2),
  - courtesy of Mark Smith and Ryan Newcombe (LHCb), with added code optimisations.

17th October 2018 John Back
* Added various decay models from LHCb EvtGenExtras package:
  - EvtBcVHad ("BC_VHAD"),
  - Evtbs2llGammaMNT ("BSTOGLLMNT"),
  - Evtbs2llGammaISRFSR ("BSTOGLLISRFSR"),
  - EvtbTosllMS ("BTOSLLMS"),
  - EvtbTosllMSExt ("BTOSLLMSEXT"),
  - EvtLb2Baryonlnu ("Lb2Baryonlnu"),
  - EvtLb2plnuLCSR ("Lb2plnuLCSR"),
  - EvtLb2plnuLQCD ("Lb2plnuLQCD"),
  - EvtFlatSqDalitz ("FLATSQDALITZ"),
  - EvtPhspFlatLifetime ("PHSPFLATLIFETIME").

5th October 2018 John Back
* Updated setupEvtGen.sh to work with the new HepForge Phabricator site.

13th March 2018 John Back
* Updated EvtPythiaEngine to correctly handle updates of various particle 
  properties so that Pythia uses the same information as EvtGen (evt.pdl)
  for the generic and alias PYTHIA decay model.

12th March 2018 John Back
* Updated EvtBcXMuNu models (X = Scalar, Vector, Tensor) to generate
  Bc to D0(star) mu nu decays, with associated form factors in EvtBCXFF, 
  - courtesy of Aleksei Luchinsky (LHCb).
* Also generalised the calculation 
  of their maximum probabilities by reusing the CalcMaxProb method in 
  EvtSemiLeptonicAmp, which now allows for different Q^2 binning 
  (default remains at 25 bins).

===
## R01-07-00

13th December 2017 John Back
* New tag incorporating all changes below.
* Recommended external packages are
  (as used in the setupEvtGen.sh script):
  - HepMC 2.06.09,
  - pythia 8.230,
  - Photos++ 3.61
  - Tauola++ 1.1.6c.

12th December 2017 John Back
* Changed Pythia validation example decay files to use Pythia8 codes.

6th December 2017 John Back
* Modified the examples to use DECAY.DEC (see 25th April 2016) instead of
  DECAY_2010.DEC. Changed EvtExternalGenList to assume Pythia8 codes are
  used in decay files by default, which is the case for DECAY.DEC. Also
  updated the setupEvtGen.sh script to work with Pythia 8.2x versions.

29th November 2017 John Back
* Modified EvtSVP, EvtVVP and EvtTVP models to handle both radiative and
  two-lepton decays,
  - courtesy of Aleksei Luchinsky (LHCb).

14th July 2017 John Back
* Only create external generator objects if they don't already exist in
  EvtExternalGenFactory.
* Modified configure script to work with Pythia 8.2x

5th July 2017 Michal Kreps
* Register the VTOSLL model.

14th June 2017 John Back
* Add isNeutralKaon() boolean function and corrected comments in EvtDDalitz.

8th May 2017 Michal Kreps
* Fix bug in EvtbTosllVectorAmp to recognise Bs --> K*bar mu mu decay as
  b --> d ll transition. 

8th May 2017 Michal Kreps
* Significantly simplify way how we decide on decay mode and daughters
  ordering in DDalitz model.
  - With new code by definition all orderings of
    daughters in the decay file will yield same output.

4th May 2017 John Back
* Further fixes to DDalitz particle ordering (including charge-conjugates):
  - Mode 5:  D0 -> K- K0bar K+ and K+ K- K0bar
  - Mode 12: D0 -> pi0 pi- pi+ and pi+ pi0 pi-
  - Removed unneeded index ordering checks for mode 10 (D+ -> pi- pi+ pi+)
    and mode 11 (Ds+ -> pi- pi+ pi+)

27th April 2017 John Back
* Fixed DDalitz particle ordering for mode 7: D+ -> pi+ K- K+ and K+ pi+ K-
  and their charge-conjugates

7th April 2017 John Back
* Modified EvtGenExternal/EvtPythiaEngine to ensure that the EvtGen-based
  instances of Pythia8 (for generic and alias decays) use the same
  particle properties as defined by EvtGen,
  - courtesy Patrick Robbe (LHCb).

5th April 2017 Michal Kreps
* Fixed indexing in copy constructor of Evt3Rank3C, which would otherwise
  produce an infinite loop;
  - bug report from David Grellscheid.

3rd November 2016 John Back
* Modified EvtFlatQ2 model to work for all B -> X lepton lepton modes, as
  well as adding an extra phase space factor to correct for the dip at low
  q^2, which is enabled by using "FLATQ2 1" instead of just "FLATQ2" in the decay file,
  - courtesy of Marcin Chrzaszcz & Thomas Blake (LHCb).


13th October 2016 John Back
* Added the TauolaCurrentOption decfile keyword to select the hadronic 
  current in Tauola; default is the BaBar-tuned current option (int = 1).
* EvtParticles can store double attributes using the functions 
  setAttributeDouble(name, double) and getAttributeDouble(name), which can
  be useful for storing and retrieving amplitude weights, for example.
  - The analogous EvtParticle integer attribute interface remains unchanged:
    setAttribute(name, int) and getAttribute(name).

13th September 2016 John Back
* Modified EvtTauolaEngine to use internal Tauola spin matrices for 
  tau pair events by temporarily setting the PDG id of the mother as a
  boson, keeping the same 4-momentum.
* The BaBar hadronic currents are now used by default.
* Also added the ability to change some Tauola parameters 
  using the "Define" keyword in decay files.
* Added an example decay file
  illustrating the new features: validation/TauolaFiles/Btautau.dec

9th September 2016 Michal Kreps
* Reimplement code in EvtBTo3pi.F, EvtBTo3piMPP.F, EvtBTo3piP00.F and 
  EvtBToKpipi.F in C++ in order to remove dependence on Fortran compiler.
  - With this, there is no internal Fortran code in EvtGen.

===
## R01-06-00

1st June 2016 John Back
* New tag incorporating all changes below.
* Recommended external packages are
  - HepMC 2.06.09
  - pythia 8.186
  - Photos++ 3.61
  - Tauola++ 1.1.5

28th April 2016 Michal Kreps
* For Ds+ --> 2pi+ pi- there was double counting of branching fraction
  resulting in total branching fraction being 1.5 times larger than measured
  one.
  - Fix by revisiting submodes, which now fill total Ds --> 3pi.

25th April 2016 Michal Kreps
* Added DECAY.DEC/XML, which contain updated semileptonic charm and beauty 
  branching fractions using the 2014 PDG, tuned to minimize disagreements 
  between measurements and EvtGen for both inclusive and exclusive decays. 
* Updated the evt.pdl particle properties file to the PDG 2014 edition.
* Implemented new LQCD form factors for Lb --> L mu mu from arXiv paper 
  1602.01399 (EvtRareLbToLllFFlQCD); old LQCD form factors are removed.

18th March 2016 John Back
* Fixed incorrect spinor algebra used in S1 -> 1/2 S2, 1/2 -> S3 S4 decays
  in EvtDiracParticle::rotateToHelicityBasis() functions,
  - courtesy of Luis Miguel Garcia Martin and the IFIC Valencia LHCb group.

19th Feburary 2016 John Back
* Fixed bug in the definition of the initial spinor term Sinit in 
  EvtRareLbToLll::HadronicAmpRS(),
  - from Tom Blake (LHCb).

12th February 2016 John Back
* From LHCb, added extensions to the EvtHQET2(FF) model for semitauonic
  decays from Brian Hamilton, which needs a patch to EvtSemiLeptonicAmp 
  from Jack Wimberley to ensure that the q^2 range is physical when 
  finding the maximum amplitude probability.

2nd December 2015 John Back
* From LHCb, added EvtKStopizmumu model for KS -> pi0 mu mu decays based on
  JHEP08(1998)004,
  - courtesy of Veronika Chobanova, Diego Martinez Santos and Jeremy Dalseno.
* Added EvtConst::Fermi for Fermi coupling constant.

===
## R01-05-00

21st October 2015 John Back
* New tag incorporating all changes below.
* Recommended external packages are
  - HepMC 2.06.09
  - pythia 8.186
  - Photos++ 3.61
  - Tauola++ 1.1.5

* Added the EvtB2MuMuMuNu model for simulating the very rare four-leptonic 
  decays B- -> mu+ mu- anti-nu_mu mu-,
  - courtesy Nikolai Nikitin.

16th October 2015 John Back
* Updated the configure script to automatically select the library names
  for PHOTOS++; version 3.56 and below uses Fortran, version 3.61 and above
  uses C++ only (default). Avoid using v3.60, since it does not work.
  This needs the PHOTOS libraries built before EvtGen is configured.
  Modified setupEvtGen.sh to use Photos++ v3.61.

7th October 2015 John Back
* Updated EvtGenExternal/EvtPhotosEngine to check that additional particles
  from the outgoing vertex are indeed (FSR) photons, since later versions of 
  PHOTOS introduce pair emission, where particles may not always be photons.
* Added the genRootDecayChain.cc validation program to create ROOT files 
  containing information about the complete decay tree. Two example test
  decay files BKstarGamma.dec and BuDst0rhop.dec can be used with this; the
  first tests PHOTOS, the second looks at sequential decay chain storage.
  The plotBKstarGamma.C ROOT macro can be used for B -> K* gamma plots.

2nd October 2015 John Back    
* Modified EvtSVPHelAmp and added a new EvtSVPHelCPMix model, implementing 
  the complete mixing phenomenology of Bs to vector gamma decays,
  - courtesy of Clara Remon (LHCb).
* EvtD0mixDalitz code: cleanup, inverted q/p for decays of D0bar (simplifies
  user decay files) and fixed y parameter bug,
  - courtesy of Jordi Tico (LHCb).
* Changed the initialisation order of the infrared cut-off in EvtPhotosEngine.
  This actually has no effect, since the exponentiation function sets it to the
  same 1e-7 value, but it is now in the correct order if we need to update it.
* Removed all remaining obsolete pragma (Win32) warnings from some classes.

23rd September 2015 Michal Kreps
* Reimplement the real Spence function in C++ and removed its fortran
  implementation.

15th September 2015 Michal Kreps
* Fixed accessed uninitialised memory in EvtPDL.cpp, line 213.
* Modified the configure and setupEvtGen.sh scripts to work on Mac; needed
  Mac compilation patch files added to the new "platform" subdirectory. 

10th September 2015 John Back
* Updated setupEvtGen.sh to use the recommended external packages:
  - HepMC 2.06.09, pythia 8.186, Photos++ 3.56 and Tauola++ 1.1.5. 
* Fixed form-factor calculations for the BTOSLLBALL model 6 used to generate
  b -> sll decays,
  - courtesy of Christoph Langenbruch and David Loh (LHCb). 
  - Affects B->K*ll, B->rholl and B->omegall, particularly the electron modes.
* In the validation directory, added runPhotosTest.sh for testing FSR in
  Upsilon(4S) -> e+ e- decays, and changed the plot comparison scripts to
  use the 2nd directory "oldRootFiles" (which could be a soft-link) for 
  including ROOT histograms made from a previous version of EvtGen.    

27th August 2015 John Back
* Added Mersenne-Twister random number generator (RNG) EvtMTRandomEngine.
  - It requires c++11 compiler features (>= gcc 4.7), which should 
    automatically be enabled by the configure script.
  - Introduced the preprocessor environment variable EVTGEN_CPP11 for c++11
    features.
  - EvtMTRandomEngine is the default RNG for the validation and test examples
    if c++11 features are enabled.
* Added a phase-space test validation/genPHSP.sh and PhaseSpacePlots.C to
  visually check the flatness of Dalitz plots in order to ensure that the
  RNG is not producing biased results that depend on particle ordering.
* Added the models EvtbsToLLLLAmp and EvtbsToLLLLHyperCP for 
  B0_q -> l+ l- l+ l- decays (SM and one supersymmetric scenario), 
  - courtesy of Nikolai Nikitin and Konstantin Toms.
  - Documentation provided in doc/evt_BQTOLLLL_model.pdf and
    doc/evt_BQTOLLLLHYPERCP_model.pdf.

* Changed the installation and set-up script name to be just setupEvtGen.sh;
  it uses the VERSION variable to specify the required tag. List of tags
  are available using either "svn ls -v http://svn.cern.ch/guest/evtgen/tags"
  or by going to http://svn.cern.ch/guest/evtgen/tags in a web browser.

12th June 2015 John Back
* Changed the width of chi_b1 in evt.pdl from 9.8928 GeV (!) to zero.

1st May 2015 John Back
* Added Bc -> scalar ell nu (EvtBcSMuNu) and Bc -> tensor ell nu 
  (EvtBcTMuNu) decays,
  - courtesy of Jack Wimberley, LHCb.
  - Also included the chi_c1 mode for EvtBcVMuNu.   

===
## R01-04-00

2nd April 2015 John Back
* Removed the EvtStdlibRandomEngine class since this can produce biases
  to kinematic distributions when one or more of the daughters is a 
  resonance, such as B0 -> K pi psi
  - (thanks to Antonio Augusto Alves Jr who discovered this issue).
  - EvtSimpleRandomEngine is now the default random number generator in the
    validation and test examples.
* Incorporated several additions and modifications from LHCb:
  a) From Michal Kreps, Tom Blake & Christoph Langenbruch,
     added rare Lb --> Lambda^(*) ell ell models described in
     arXiv:1108.6129, with various form factors from Gutsche et al. 
     (arXiv:1301.3737) and lattice QCD  (arXiv:1212.4827)
  b) From Andrew Crocombe, implemented Bs --> K* form factors 
     from Ball-Zwicky and z-parametrization form factors from 
     arXiv:1006.4945 for EvtbTosllBallFF
  c) Christoph Langenbruch fixed the Bs -> phi ll form factors in
     EvtbTosllBallFF; T3 showed a non-physical pole at very low 
     q2 which significantly affected the electron mode
  d) From Michal Kreps, removed semicolons from wrong places to
     clear warnings when compiled with the "-pedantic" option.

9th October 2014 John Back
* Change svnweb.cern.ch to svn.cern.ch in the setup script.

1st April 2014 John Back
* In EvtReport, modified the logging output severity status flags
  to have the "EVTGEN_" prefix, e.g. INFO becomes EVTGEN_INFO. 
* The global report() function has been renamed to EvtGenReport().

31st March 2014 John Back
* Added the ability to store named attributes for EvtParticles in the
  form of a map<string, int>. The setAttribute(name, value) stores the
  required value, while getAttribute(name) retrieves the integer value.
  This is used in EvtPhotosEngine to specify the final-state radiation
  "FSR" attribute to 1 for any additional photons (EvtPhotonParticles)
  created by Photos++. It also stores the "ISR" attribute, but this 
  is always set to zero, since only FSR photons are created.
  If the named attribute does not exist, then getAttribute() returns zero.

29th January 2014 Daniel Craik
* Removed mass assertion on GS shape in EvtDalitzReso to allow it to also 
  be used for charged rho resonances.

27th January 2014 John Back
* Minor corrections to Vub models to remove further gcc 4.8 warnings.
* Updated configure script to work for MacOS clang (from Genser team).

===
## R01-03-00

9th January 2014 John Back
* New tag version "1.3.0", incorporating all changes below.
* Replaced auto-install script to work with this version as well as
  the latest versions of all external generator packages.
* Updated README to mention the new CERN-based web pages for Photos++ 
  and Tauola++.

8th January 2014 John Back
* Fix gcc 4.6 and 4.8 compilation warnings,
  - courtesy of Patrick Robbe (LHCb);
  - main changes are removal of unused variables.
* Changed the EvtPythiaEngine class and configure script to use new 
  Pythia 8 header locations; Pythia 8.180 or above is now required.

7th January 2014 John Back
* Modified EvtBCVFF to correct the Kiselev form factors
  - from Jack Wimberley (LHCb).

9th October 2013 Daniel Craik
* Added Gounaris-Sakurai and Gaussian shapes to EvtGenericDalitz 
  and set sensible defaults for LASS parameters.

19th September 2013 John Back
* Modified EvtGenExternal/EvtPythiaEngine to keep track of any new 
  particles that are added to the default Pythia database to avoid
  duplicating particle/anti-particle entries that could override
  previously defined Pythia decay chains.

18th September 2013 John Back
* Added Mac OS flags for the configure script and src/Makefile.

15th July 2013 Daniel Craik
* Added flag to turn on scaling of LASS amplitude by M/q in EvtDalitzReso

15th July 2013 Daniel Craik
* EvtParserXML now accepts file names containing environment variables,
  exponential non-resonant shape in EvtDalitzReso now defined as exp(-alpha*m^2),
  LASS shape in EvtDalitzReso now takes a cutoff parameter

4th July 2013 Daniel Craik
* Added LASS, exponential non-resonant and linear non-resonant shapes to EvtGenericDalitz.

3rd July 2013 Daniel Craik
* Fixed auto-install script for R01-02-00.

1st July 2013 Daniel Craik
* Added auto-install script for R01-02-00.

===
## R01-02-00

15th May 2013 John Back
* New tag, version "1.2.0", incorporating all changes below.

14th May 2013 Michal Kreps
* Added Blatt-Weisskopf barrier factors up to L=5 in 
  EvtGenBase/EvtBlattWeisskopf::compute().

14th May 2013 John Back
* Added additional entries (appended at the end) to the evt.pdl particle 
  data file
  - courtesy of Romulus Godang and Belle II colleagues.

14th March 2013 John Back
* Added the method EvtParticle::getPDGId() to get the PDG integer for a
  particle directly (which just calls EvtPDL::getStdHep()).
* Added a check in EvtPhotosEngine::doDecay to skip Photos if a given 
  particle has too many daughters (>= 10) to avoid a problem with a 
  hard coded upper limit in the Photos PHOENE subroutine.

2nd February 2013 Daniel Craik
* Updated EvtDalitzTable to estimate probMax when it is missing from a 
  Dalitz model.

1st  February 2013 John Back
* Added the ability to read in Pythia 6 commands in ascii decay files in
  EvtDecayTable::readDecayFile (this was already possible in xml files).
* Modified the Photos++ engine default settings to be more suited to B
  decays (from LHCb defaults).

31st January 2013 John Back
* Added the ability to read in Pythia 8 commands in ascii decay files
  in EvtDecayTable::readDecayFile. They can be set using the syntax:
  "PythiaTypeParam module:variable=value", where Type = Generic, Alias or
  Both for specifying whether the parameter is for the generic or alias
  Pythia decay engines (or both). The 2nd argument must not contain any
  spaces. Fixed the list of commands strings being used in the 
  EvtPythiaEngine class (i.e. Pythia parameters that can be set via decay
  files).

31st January 2013 Daniel Craik
* Added named parameters to various decay models.

30th January 2013 John Back
* Fixed some of the parameter arguments used in the EvtSVSCPiso model.

24th January 2013 John Back
* Set the Photos++ and Tauola++ models to use the EvtGen random number 
  engine when useEvtGenRandom is set to true in the EvtExternalGenList 
  constructor.

23rd January 2013 John Back
* Added EvtGenExternal/EvtPythiaRandom to allow the use of the EvtGen
  random number engine to also be used for the random engine for Pythia 8.
* Added a boolean (useEvtGenRandom, default = true) within the 
  EvtExternalGenList constructor to use this feature.

18th December 2012 John Back
* Corrected some wrong daughter ordering assignments for decay modes 5 and
  12 in EvtDDalitz. Updated validation/DalitzDecays.xml to also contain
  D decay mode 12, as well as various modes that may use K_S0 and K_L0. 
* Added validation/genDDalitzModes.sh and updated validation/compareDalitz.C to
  do a complete comparison of the D Dalitz modes with the xml versions.

11th December 2012 Daniel Craik
* Updated the Xml parser to support named model parameters.
* Updated the generic Dalitz model to use named model parameters as an example.

15th October 2012 John Back
* Make EvtSimpleRandomEngine inherit from EvtRandomEngine to avoid
  crash in EvtGen.cpp when no random engine is defined 
  - (from Bjoern Spruck).

===
## R01-01-00

4th October 2012 John Back
* New tag, version "1.1.0", incorporating all changes below.

* Provide proper default constructors for EvtVector4R and 
  EvtPhotonParticle. Modified the validation and test code to also
  compile/link in the case of no external generators being included.

3rd October 2012 John Back
* Corrected the t3 vector form factor values for the Ball-Zwicky 2005 
  model (modelId = 6) in EvtbTosllBallFF::getVectorFF(), which
  were set to t3tilde instead.

18th September 2012 John Back
* Moved the external generator engines to a new sub-directory 
  EvtGenExternal. Building the code now creates 2 libraries:
  libEvtGen.so (Base+Models) and libEvtGenExternal.so.
  - This allows anyone to ignore using the new external generators
    if required (by not creating/loading the 2nd library).
* Added prefix option to the configure script/Makefile to allow the user 
  to specify an installation directory for the include files, libraries, 
  DECAY.DEC and evt.pdl files (for Genser).

14th September 2012 Michal Kreps
* Fixed the calculation of the angle between decay planes in the function
  EvtKine::EvtDecayAngleChi. Fixed typo in EvtLb2Lll decay model. Only 
  some NP scenarious could be affected, SM one is definitely unaffected.

13th September 2012 John Back
* Added the use of the environment variables EVTGEN_PHOTOS, EVTGEN_PYTHIA 
  and EVTGEN_TAUOLA to specify if the Photos, Pythia and/or Tauola engine
  classes are used or not. These variables are set by the configure script,
  depending if the library paths are specified for these generators.

===
## R01-00-01

12th September 2012 John Back
* New tag incorporating all changes below, since R01-00-00.

11th September 2012 John Back
* Modified the Photos and Tauola engine classes to use the new
  Photospp and Tauolapp namespaces that are present in the latest
  versions of Photos++(3.5) and Tauola++(1.0.7).
* Updated the configure file to get the correct location of the Tauola++
  include files. 
* Added the D0->pi+pi-pi0 decay mode in EvtDDalitz
  - from Marco Gersabeck and Frederic Dreyer (LHCb).
* Added new decay models/classes from Alexey Luchinsky (LHCb):
  EvtBcVMuNu, EvtTVP, EvtWnPi, EvtSVP, EvtXPsiGamma, EvtBcVNpi

29th June 2012 John Back
* Corrected mass(squared) variables filled in the Dalitz TTree in 
  validation/genExampleRootFiles.

15th May 2012 Daniel Craik
* Updated EvtD0gammaDalitz to deal with D mesons from neutral B->DK
* Added save function to validation/compareDalitz.C.

11th May 2012 Daniel Craik
* Replaced BaBar specific configuration for BlattWeisskopf birth factors.
* Updated XML conversion script to handle new configuration.
* Fixed some bugs in the XML conversion script related to particle 
  modifications.

9th May 2012 Daniel Craik
* Added latex documentation for xml decay files.

2nd May 2012 Daniel Craik
* Added resDaughters attribute to the Dalitz resonance xml tag to 
  simplify defining symmetric resonances. Updated validation xml files to 
  use the new functionality.

27th April 2012 Daniel Craik
* Upgraded EvtGenericDalitz to use EvtDalitzReso for resonances.
* Added validation to compare EvtGenericDalitz to all 11 EvtDDalitz modes.
* Added a root macro to quickly compare two Dalitz decays for validation.

24th April 2012 John Back
* Solved two bugs in the EvtD0gammaDalitz model (from Jordi Tico, LHCb):
  configuration of the conjugated model, and using only the B charge
  to determine the model used, not the D flavour.

17th April 2012 Daniel Craik
* Updated the GenericDalitz validation code to use the same probMax 
  values as DDalitz.
* Added XML decay file parsing to EvtGen::readUDecay.
  - Dec files are still the default.

30th March 2012 John Back
* Update maximum probability values in EvtDDalitz::initProbMax()
  for all DDalitz modes.

23rd March 2012 John Back
* Added the EvtEta2MuMuGamma decay model from LHCb.

21st March 2012 John Back
* Added EvtD0gammaDalitz decay model from LHCb.

20th March 2012 Daniel Craik
* Added backwards compatibility for Pythia 6 commands in the XML configuration.
* Updated decay file conversion tool to convert JetSetPar lines to pythia6Param
  tags.

19th March 2012 Daniel Craik
* Added infrastructure to pass commands to external generators.
* XML config now takes Pythia8 configuration commands.

16th March 2012 Daniel Craik
* Added the ability to define particles from the PDL for Dalitz decay
  resonances instead of defining mass, width and spin seperately.
* Renamed the lifetime attribute of Dalitz decay resonaces to width to avoid
  confusion.
* Added further validation code for the generic Dalitz model.

15th March 2012 Daniel Craik
* Added validation code for xml decay files and the generic Dalitz model. 

===
## R01-00-00 

6th March 2012 John Back
* First official version for Genser (evtgen 1.0.0) that includes
 support for external generators: Pythia8, Photos++ and Tauola++. 
* This also includes a preliminary version of creating Dalitz plot
 decay models using EvtGenericDalitz.
