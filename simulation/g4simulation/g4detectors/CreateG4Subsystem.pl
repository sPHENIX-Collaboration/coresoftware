#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

sub CreateSubsystemInclude;
sub CreateSubsystemImplementation;

sub CreateDetectorInclude;
sub CreateDetectorImplementation;

sub CreateSteppingActionInclude;
sub CreateSteppingActionImplementation;

sub CreateMakefile;
sub CreateAutogen;
sub CreateConfigure;

if ($#ARGV < 0)
{
    print "Repo + Location: coresoftware:simulation/g4simulation/g4detectors\n";
    print "Usage:\n";
    print "CreateG4Subsystem.pl <Detector Name>\n";
    print "options:\n";
    print "--all : create also autogen.sh, configure.ac and Makefile.am\n";
    print "--overwrite : overwrite existing files\n";
    exit(0);
}

my $createall;
my $overwrite;
GetOptions("all" => \$createall, "overwrite" => \$overwrite);
my $detectorname = sprintf("%s",$ARGV[0]);
if ($detectorname =~ m/[^a-zA-Z0-9]/)
{
    print "Detector name contains invalid characters - allowed are alphanumeric (lower and upper caps)\n";
    exit(1);
}
if (substr($detectorname,0,1) =~ m/[0-9]/)
{
    print "Detector name must start with a letter\n";
    exit(1);
}

my $subsysclassname = sprintf("%sSubsystem",$detectorname);
my $detectorclassname = sprintf("%sDetector",$detectorname);
my $steppingclassname = sprintf("%sSteppingAction",$detectorname);

my %listoffiles = ();

my $subsystem_includefile = sprintf("%s.h",$subsysclassname);
my $subsystem_implementfile = sprintf("%s.cc",$subsysclassname );
my $detector_includefile = sprintf("%s.h", $detectorclassname);
my $detector_implementfile = sprintf("%s.cc", $detectorclassname);
my $steppingaction_includefile = sprintf("%s.h", $steppingclassname);
my $steppingaction_implementfile = sprintf("%s.cc", $steppingclassname);

$listoffiles{$subsystem_includefile} = $detectorname;
$listoffiles{$subsystem_implementfile} = $detectorname;
$listoffiles{$detector_includefile} = $detectorname;
$listoffiles{$detector_implementfile} = $detectorname;
$listoffiles{$steppingaction_includefile} = $detectorname;
$listoffiles{$steppingaction_implementfile} = $detectorname;

if (defined $createall)
{
    $listoffiles{"autogen.sh"} = $detectorname;
    $listoffiles{"configure.ac"} = $detectorname;
    $listoffiles{"Makefile.am"} = $detectorname;
}

# check if files exist if overwrite is not set

if (! defined $overwrite)
{
    foreach my $file (keys %listoffiles)
    {
	if (-f $file)
	{
	    print "$file exists but overwrite option not set\n";
	    exit(1);
	}

    }
}

CreateSubsystemInclude($subsystem_includefile);
CreateSubsystemImplementation($subsystem_implementfile);

CreateDetectorInclude($detector_includefile);
CreateDetectorImplementation($detector_implementfile);

CreateSteppingActionInclude($steppingaction_includefile);
CreateSteppingActionImplementation($steppingaction_implementfile);

if (defined $createall)
{
    CreateAutogen();
    CreateMakefile();
    CreateConfigure();
}
exit(0);

sub CreateSteppingActionInclude()
{
    my $file = shift;
    open(F,">$file");
    my $includeguard = uc(sprintf("%s_H",$steppingclassname));
    open(F,">$file");
    print F "// Tell emacs that this is a C++ source\n";
    print F "//  -*- C++ -*-.\n";
    print F "#ifndef $includeguard\n";
    print F "#define $includeguard\n";
    print F "\n";
    print F "#include <g4main/PHG4SteppingAction.h>\n";
    print F "\n";

    print F "class $detectorclassname;\n";
    print F "\n";

    print F "class G4Step;\n";
    print F "class G4VPhysicalVolume;\n";
    print F "class PHCompositeNode;\n";
    print F "class PHG4Hit;\n";
    print F "class PHG4HitContainer;\n";
    print F "class PHParameters;\n";
    print F "\n";

    print F "class $steppingclassname : public PHG4SteppingAction\n";
    print F "{\n";
    print F " public:\n";
    print F "  //! constructor\n";
    print F "  $steppingclassname($detectorclassname*, const PHParameters* parameters);\n";
    print F "\n";

    print F "  //! destructor\n";
    print F "  virtual ~$steppingclassname();\n";
    print F "\n";

    print F "  //! stepping action\n";
    print F "  virtual bool UserSteppingAction(const G4Step*, bool);\n";
    print F "\n";

    print F "  //! reimplemented from base class\n";
    print F "  virtual void SetInterfacePointers(PHCompositeNode*);\n";
    print F "\n";

    print F " private:\n";
    print F "  //! pointer to the detector\n";
    print F "  $detectorclassname* m_Detector;\n";
    print F "  const PHParameters* m_Params;\n";
    print F "  //! pointer to hit container\n";
    print F "  PHG4HitContainer* m_HitContainer;\n";
    print F "  PHG4Hit* m_Hit;\n";
    print F "  PHG4HitContainer* m_SaveHitContainer;\n";
    print F "  G4VPhysicalVolume* m_SaveVolPre;\n";
    print F "  G4VPhysicalVolume* m_SaveVolPost;\n";
    print F "\n";

    print F "  int m_SaveTrackId;\n";
    print F "  int m_SavePreStepStatus;\n";
    print F "  int m_SavePostStepStatus;\n";
    print F "  int m_ActiveFlag;\n";
    print F "  int m_BlackHoleFlag;\n";
    print F "  double m_EdepSum;\n";
    print F "  double m_EionSum;\n";
    print F "};\n";
    print F "\n";

    print F "#endif // $includeguard\n";
    close(F);
}

sub CreateSteppingActionImplementation()
{
    my $file = shift;
    open(F,">$file");
    print F "//____________________________________________________________________________..\n";
    print F "//\n";
    print F "// This is a working template for the Stepping Action which needs to be implemented\n";
    print F "// for active detectors. Most of the code is error handling and access to the G4 objects\n";
    print F "// and our data structures. It does not need any adjustment. The only thing you need to\n";
    print F "// do is to add the properties of the G4Hits you want to save for later analysis\n";
    print F "// This needs to be done in 2 places, G4Hits are generated when a G4 track enters a new\n";
    print F "// volume (or is created). Here you give it an initial value. When the G4 track leaves\n";
    print F "// the volume the final value needs to be set.\n";
    print F "// The places to do this is marked by //implement your own here//\n";
    print F "//\n";

    print F "// As guidance you can look at the total (integrated over all steps in a volume) energy\n";
    print F "// deposit which should always be saved.\n";
    print F "// Additionally the total ionization energy is saved - this can be removed if you are not\n";
    print F "// interested in this. Naturally you may want remove these comments in your version\n";
    print F "//\n";
    print F "//____________________________________________________________________________..\n";
    print F "\n";
    print F "#include \"$steppingaction_includefile\"\n";
    print F "\n";

    print F "#include \"$detector_includefile\"\n";
    print F "\n";

    print F "#include <phparameter/PHParameters.h>\n";
    print F "\n";

    print F "#include <g4detectors/PHG4StepStatusDecode.h>\n";
    print F "\n";

    print F "#include <g4main/PHG4Hit.h>\n";
    print F "#include <g4main/PHG4HitContainer.h>\n";
    print F "#include <g4main/PHG4Hitv1.h>\n";
    print F "#include <g4main/PHG4Shower.h>\n";
    print F "#include <g4main/PHG4SteppingAction.h>\n";
    print F "#include <g4main/PHG4TrackUserInfoV1.h>\n";
    print F "\n";

    print F "#include <phool/getClass.h>\n";
    print F "\n";

    print F "#include <TSystem.h>\n";
    print F "\n";

    print F "#include <Geant4/G4ParticleDefinition.hh> \n";
    print F "#include <Geant4/G4ReferenceCountedHandle.hh>\n";
    print F "#include <Geant4/G4Step.hh>\n";
    print F "#include <Geant4/G4StepPoint.hh> \n";
    print F "#include <Geant4/G4StepStatus.hh>\n";
    print F "#include <Geant4/G4String.hh> \n";
    print F "#include <Geant4/G4SystemOfUnits.hh>\n";
    print F "#include <Geant4/G4ThreeVector.hh>\n";
    print F "#include <Geant4/G4TouchableHandle.hh>\n";
    print F "#include <Geant4/G4Track.hh>\n";
    print F "#include <Geant4/G4TrackStatus.hh>\n";
    print F "#include <Geant4/G4Types.hh>\n";
    print F "#include <Geant4/G4VPhysicalVolume.hh>\n";
    print F "#include <Geant4/G4VTouchable.hh>\n";
    print F "#include <Geant4/G4VUserTrackInformation.hh>\n";
    print F "\n";

    print F "#include <cmath>\n";
    print F "#include <iostream>\n";
    print F "#include <string>\n";
    print F "\n";

    print F "class PHCompositeNode;\n";
    print F "\n";

    print F "//____________________________________________________________________________..\n";
    print F "$steppingclassname\:\:$steppingclassname($detectorclassname *detector, const PHParameters *parameters)\n";
    print F "  : PHG4SteppingAction(detector->GetName())\n";
    print F "  , m_Detector(detector)\n";
    print F "  , m_Params(parameters)\n";
    print F "  , m_HitContainer(nullptr)\n";
    print F "  , m_Hit(nullptr)\n";
    print F "  , m_SaveHitContainer(nullptr)\n";
    print F "  , m_SaveVolPre(nullptr)\n";
    print F "  , m_SaveVolPost(nullptr)\n";
    print F "  , m_SaveTrackId(-1)\n";
    print F "  , m_SavePreStepStatus(-1)\n";
    print F "  , m_SavePostStepStatus(-1)\n";
    print F "  , m_ActiveFlag(m_Params->get_int_param(\"active\"))\n";
    print F "  , m_BlackHoleFlag(m_Params->get_int_param(\"blackhole\"))\n";
    print F "  , m_EdepSum(0)\n";
    print F "  , m_EionSum(0)\n";
    print F "{\n";
    print F "}\n";
    print F "\n";

    print F "//____________________________________________________________________________..\n";
    print F "$steppingclassname\:\:~$steppingclassname()\n";
    print F "{\n";
    print F "  // if the last hit was a zero energie deposit hit, it is just reset\n";
    print F "  // and the memory is still allocated, so we need to delete it here\n";
    print F "  // if the last hit was saved, hit is a nullptr pointer which are\n";
    print F "  // legal to delete (it results in a no operation)\n";
    print F "  delete m_Hit;\n";
    print F "}\n";
    print F "\n";

    print F "//____________________________________________________________________________..\n";
    print F "// This is the implementation of the G4 UserSteppingAction\n";
    print F "bool $steppingclassname\:\:UserSteppingAction(const G4Step *aStep,bool was_used)\n";
    print F "{\n";
    print F "  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();\n";
    print F "  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();\n";
    print F "  // get volume of the current step\n";
    print F "  G4VPhysicalVolume *volume = touch->GetVolume();\n";
    print F "  // IsInDetector(volume) returns\n";
    print F "  //  == 0 outside of detector\n";
    print F "  //   > 0 for hits in active volume\n";
    print F "  //  < 0 for hits in passive material\n";
    print F "  int whichactive = m_Detector->IsInDetector(volume);\n";
    print F "  if (!whichactive)\n";
    print F "  {\n";
    print F "    return false;\n";
    print F "  }\n";

    print F "  // collect energy and track length step by step\n";
    print F "  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;\n";
    print F "  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;\n";
    print F "  const G4Track *aTrack = aStep->GetTrack();\n";

    print F "  // if this detector stops everything, just put all kinetic energy into edep\n";
    print F "  if (m_BlackHoleFlag)\n";
    print F "  {\n";
    print F "    edep = aTrack->GetKineticEnergy() / GeV;\n";
    print F "    G4Track *killtrack = const_cast<G4Track *>(aTrack);\n";
    print F "    killtrack->SetTrackStatus(fStopAndKill);\n";
    print F "  }\n";

    print F "  // we use here only one detector in this simple example\n";
    print F "  // if you deal with multiple detectors in this stepping action\n";
    print F "  // the detector id can be used to distinguish between them\n";
    print F "  // hits can easily be analyzed later according to their detector id\n";
    print F "  int detector_id = 0;  // we use here only one detector in this simple example\n";
    print F "  bool geantino = false;\n";
    print F "  // the check for the pdg code speeds things up, I do not want to make\n";
    print F "  // an expensive string compare for every track when we know\n";
    print F "  // geantino or chargedgeantino has pid=0\n";
    print F "  if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&\n";
    print F "      aTrack->GetParticleDefinition()->GetParticleName().find(\"geantino\") !=\n";
    print F "          std::string::npos)  // this also accounts for \"chargedgeantino\"\n";
    print F "  {\n";
    print F "    geantino = true;\n";
    print F "  }\n";
    print F "  G4StepPoint *prePoint = aStep->GetPreStepPoint();\n";
    print F "  G4StepPoint *postPoint = aStep->GetPostStepPoint();\n";
    print F "\n";
    print F "// Here we have to decide if we need to create a new hit.  Normally this should \n";
    print F "// only be neccessary if a G4 Track enters a new volume or is freshly created\n";
    print F "// For this we look at the step status of the prePoint (beginning of the G4 Step).\n";
    print F "// This should be either fGeomBoundary (G4 Track crosses into volume) or \n";
    print F "// fUndefined (G4 Track newly created)\n";
    print F "// Sadly over the years with different G4 versions we have observed cases where\n";
    print F "// G4 produces \"impossible hits\" which we try to catch here\n";
    print F "// These errors were always rare and it is not clear if they still exist but we\n";
    print F "// still check for them for safety. We can reproduce G4 runs identically (if given\n";
    print F "// the sequence of random number seeds you find in the log), the printouts help\n";
    print F "// us giving the G4 support information about those failures\n";
    print F "// \n";
    print F "  switch (prePoint->GetStepStatus())\n";
    print F "  {\n";
    print F "  case fPostStepDoItProc:\n";
    print F "    if (m_SavePostStepStatus != fGeomBoundary)\n";
    print F "    {\n";
    print F "      // this is the okay case, fPostStepDoItProc called in a volume, not first thing inside\n";
    print F "      // a new volume, just proceed here\n";
    print F "      break;\n";
    print F "    }\n";
    print F "    else\n";
    print F "    {\n";
    print F "      // this is an impossible G4 Step print out diagnostic to help debug, not sure if\n";
    print F "      // this is still with us\n";
    print F "      std::cout << GetName() << \": New Hit for  \" << std::endl;\n";
    print F "      std::cout << \"prestep status: \"\n";
    print F "           << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())\n";
    print F "           << \", poststep status: \"\n";
    print F "           << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())\n";
    print F "           << \", last pre step status: \"\n";
    print F "           << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)\n";
    print F "           << \", last post step status: \"\n";
    print F "           << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;\n";
    print F "      std::cout << \"last track: \" << m_SaveTrackId\n";
    print F "           << \", current trackid: \" << aTrack->GetTrackID() << std::endl;\n";
    print F "      std::cout << \"phys pre vol: \" << volume->GetName()\n";
    print F "           << \" post vol : \" << touchpost->GetVolume()->GetName() << std::endl;\n";
    print F "      std::cout << \" previous phys pre vol: \" << m_SaveVolPre->GetName()\n";
    print F "           << \" previous phys post vol: \" << m_SaveVolPost->GetName() << std::endl;\n";
    print F "    }\n";
    print F "// These are the normal cases\n";
    print F "  case fGeomBoundary:\n";
    print F "  case fUndefined:\n";
    print F "    if (!m_Hit)\n";
    print F "    {\n";
    print F "      m_Hit = new PHG4Hitv1();\n";
    print F "    }\n";
    print F "    m_Hit->set_layer(detector_id);\n";
    print F "    // here we set the entrance values in cm\n";
    print F "    m_Hit->set_x(0, prePoint->GetPosition().x() / cm);\n";
    print F "    m_Hit->set_y(0, prePoint->GetPosition().y() / cm);\n";
    print F "    m_Hit->set_z(0, prePoint->GetPosition().z() / cm);\n";
    print F "    // time in ns\n";
    print F "    m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);\n";
    print F "    // set the track ID\n";
    print F "    m_Hit->set_trkid(aTrack->GetTrackID());\n";
    print F "    m_SaveTrackId = aTrack->GetTrackID();\n";
    print F "    // set the initial energy deposit\n";
    print F "    m_EdepSum = 0;\n";
    print F "    // implement your own here://\n";
    print F "    // add the properties you are interested in via set_XXX methods\n";
    print F "    // you can find existing set methods in \$OFFLINE_MAIN/include/g4main/PHG4Hit.h\n";
    print F "    // this is initialization of your value. This is not needed you can just set the final\n";
    print F "    // value at the last step in this volume later one\n";
    print F "    if (whichactive > 0)\n";
    print F "    {\n";
    print F "      m_EionSum = 0;  // assuming the ionization energy is only needed for active\n";
    print F "                      // volumes (scintillators)\n";
    print F "      m_Hit->set_eion(0);\n";
    print F "      m_SaveHitContainer = m_HitContainer;\n";
    print F "    }\n";
    print F "    else\n";
    print F "    {\n";
    print F "      std::cout << \"implement stuff for whichactive < 0 (inactive volumes)\" << std::endl;\n";
    print F "      gSystem->Exit(1);\n";
    print F "    }\n";
    print F "    // this is for the tracking of the truth info\n";
    print F "    if (G4VUserTrackInformation *p = aTrack->GetUserInformation())\n";
    print F "    {\n";
    print F "      if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))\n";
    print F "      {\n";
    print F "        m_Hit->set_trkid(pp->GetUserTrackId());\n";
    print F "        pp->GetShower()->add_g4hit_id(m_SaveHitContainer->GetID(), m_Hit->get_hit_id());\n";
    print F "      }\n";
    print F "    }\n";

    print F "    break;\n";
    print F "  default:\n";
    print F "    break;\n";
    print F "  }\n";
    print F "\n";

    print F "  // This section is called for every step\n";
    print F "  // some sanity checks for inconsistencies (aka bugs) we have seen over the years\n";
    print F "  // check if this hit was created, if not print out last post step status\n";
    print F "  if (!m_Hit || !std::isfinite(m_Hit->get_x(0)))\n";
    print F "  {\n";
    print F "    std::cout << GetName() << \": hit was not created\" << std::endl;\n";
    print F "    std::cout << \"prestep status: \"\n";
    print F "         << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())\n";
    print F "         << \", poststep status: \"\n";
    print F "         << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())\n";
    print F "         << \", last pre step status: \"\n";
    print F "         << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)\n";
    print F "         << \", last post step status: \"\n";
    print F "         << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;\n";
    print F "    std::cout << \"last track: \" << m_SaveTrackId\n";
    print F "         << \", current trackid: \" << aTrack->GetTrackID() << std::endl;\n";
    print F "    std::cout << \"phys pre vol: \" << volume->GetName()\n";
    print F "         << \" post vol : \" << touchpost->GetVolume()->GetName() << std::endl;\n";
    print F "    std::cout << \" previous phys pre vol: \" << m_SaveVolPre->GetName()\n";
    print F "         << \" previous phys post vol: \" << m_SaveVolPost->GetName() << std::endl;\n";
    print F "    // This is fatal - a hit from nowhere. This needs to be looked at and fixed\n";
    print F "    gSystem->Exit(1);\n";
    print F "  }\n";
    print F "  // check if track id matches the initial one when the hit was created\n";
    print F "  if (aTrack->GetTrackID() != m_SaveTrackId)\n";
    print F "  {\n";
    print F "    std::cout << GetName() << \": hits do not belong to the same track\" << std::endl;\n";
    print F "    std::cout << \"saved track: \" << m_SaveTrackId\n";
    print F "         << \", current trackid: \" << aTrack->GetTrackID()\n";
    print F "         << \", prestep status: \" << prePoint->GetStepStatus()\n";
    print F "         << \", previous post step status: \" << m_SavePostStepStatus << std::endl;\n";
    print F "    // This is fatal - a hit from nowhere. This needs to be looked at and fixed\n";
    print F "    gSystem->Exit(1);\n";
    print F "  }\n";
    print F "\n";

    print F "// We need to cache a few things from one step to the next\n";
    print F "// to identify impossible hits and subsequent debugging printout\n";
    print F "  m_SavePreStepStatus = prePoint->GetStepStatus();\n";
    print F "  m_SavePostStepStatus = postPoint->GetStepStatus();\n";
    print F "  m_SaveVolPre = volume;\n";
    print F "  m_SaveVolPost = touchpost->GetVolume();\n";

    print F "  // here we just update the exit values, it will be overwritten\n";
    print F "  // for every step until we leave the volume or the particle\n";
    print F "  // ceases to exist\n";
    print F "  // sum up the energy to get total deposited\n";
    print F "  m_EdepSum += edep;\n";
    print F "  if (whichactive > 0)\n";
    print F "  {\n";
    print F "    m_EionSum += eion;\n";
    print F "  }\n";
    print F "  // if any of these conditions is true this is the last step in\n";
    print F "  // this volume and we need to save the hit\n";
    print F "  // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume\n";
    print F "  // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world\n";
    print F "  // (happens when your detector goes outside world volume)\n";
    print F "  // postPoint->GetStepStatus() == fAtRestDoItProc: track stops (typically\n";
    print F "  // aTrack->GetTrackStatus() == fStopAndKill is also set)\n";
    print F "  // aTrack->GetTrackStatus() == fStopAndKill: track ends\n";
    print F "  if (postPoint->GetStepStatus() == fGeomBoundary ||\n";
    print F "      postPoint->GetStepStatus() == fWorldBoundary ||\n";
    print F "      postPoint->GetStepStatus() == fAtRestDoItProc ||\n";
    print F "      aTrack->GetTrackStatus() == fStopAndKill)\n";
    print F "  {\n";
    print F "    // save only hits with energy deposit (or geantino)\n";
    print F "    if (m_EdepSum > 0 || geantino)\n";
    print F "    {\n";
    print F "      // update values at exit coordinates and set keep flag\n";
    print F "      // of track to keep\n";
    print F "      m_Hit->set_x(1, postPoint->GetPosition().x() / cm);\n";
    print F "      m_Hit->set_y(1, postPoint->GetPosition().y() / cm);\n";
    print F "      m_Hit->set_z(1, postPoint->GetPosition().z() / cm);\n";

    print F "      m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);\n";
    print F "      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())\n";
    print F "      {\n";
    print F "        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))\n";
    print F "        {\n";
    print F "          pp->SetKeep(1);  // we want to keep the track\n";
    print F "        }\n";
    print F "      }\n";
    print F "      if (geantino)\n";
    print F "      {\n";
    print F " //implement your own here://\n";
    print F " // if you want to do something special for geantinos (normally you do not)\n";
    print F "        m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way\n";
    print F "                              // geantinos survive the g4hit compression\n";
    print F "        if (whichactive > 0)\n";
    print F "        {\n";
    print F "          m_Hit->set_eion(-1);\n";
    print F "        }\n";
    print F "      }\n";
    print F "      else\n";
    print F "      {\n";
    print F "        m_Hit->set_edep(m_EdepSum);\n";
    print F "      }\n";
    print F " //implement your own here://\n";
    print F " // what you set here will be saved in the output\n";
    print F "      if (whichactive > 0)\n";
    print F "      {\n";
    print F "        m_Hit->set_eion(m_EionSum);\n";
    print F "      }\n";
    print F "      m_SaveHitContainer->AddHit(detector_id, m_Hit);\n";
    print F "      // ownership has been transferred to container, set to null\n";
    print F "      // so we will create a new hit for the next track\n";
    print F "      m_Hit = nullptr;\n";
    print F "    }\n";
    print F "    else\n";
    print F "    {\n";
    print F "      // if this hit has no energy deposit, just reset it for reuse\n";
    print F "      // this means we have to delete it in the dtor. If this was\n";
    print F "      // the last hit we processed the memory is still allocated\n";
    print F "      m_Hit->Reset();\n";
    print F "    }\n";
    print F "  }\n";
    print F "  // return true to indicate the hit was used\n";
    print F "  return true;\n";
    print F "}\n";
    print F "\n";

    print F "//____________________________________________________________________________..\n";
    print F "void $steppingclassname\:\:SetInterfacePointers(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  std::string hitnodename = \"G4HIT_\" + m_Detector->GetName();\n";

    print F "  // now look for the map and grab a pointer to it.\n";
    print F "  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);\n";

    print F "  // if we do not find the node we need to make it.\n";
    print F "  if (!m_HitContainer)\n";
    print F "  {\n";
    print F "    std::cout << \"$steppingclassname\:\:SetTopNode - unable to find \"\n";
    print F "              << hitnodename << std::endl;\n";
    print F "  }\n";
    print F "}\n";
    close(F);
}

sub CreateDetectorImplementation()
{
    my $file = shift;
    open(F,">$file");
    print F "//____________________________________________________________________________..\n";
    print F "//\n";
    print F "// This is a working template for the G4 Construct() method which needs to be implemented\n";
    print F "// We wedge a method between the G4 Construct() to enable volume hierarchies on the macro\n";
    print F "// so here it is called ConstructMe() but there is no functional difference\n";
    print F "// Currently this installs a simple G4Box solid, creates a logical volume from it\n";
    print F "// and places it. Put your own detector in place (just make sure all active volumes\n";
    print F "// get inserted into the m_PhysicalVolumesSet)\n";
    print F "// \n";
    print F "// Rather than using hardcoded values you should consider using the parameter class\n";
    print F "// Parameter names and defaults are set in $subsysclassname\:\:SetDefaultParameters()\n";
    print F "// Only parameters defined there can be used (also to override in the macro)\n";
    print F "// to avoids typos.\n";
    print F "// IMPORTANT: parameters have no inherent units, there is a convention (cm/deg)\n";
    print F "// but in any case you need to multiply them here with the correct CLHEP/G4 unit \n";
    print F "// \n";
    print F "// The place where you put your own detector is marked with\n";
    print F "// //begin implement your own here://\n";
    print F "// //end implement your own here://\n";
    print F "// Do not forget to include the G4 includes for your volumes\n";
    print F "//____________________________________________________________________________..\n";
    print F "\n";
    print F "#include \"$detector_includefile\"\n";
    print F "\n";
    print F "#include <phparameter/PHParameters.h>\n";
    print F "\n";
    print F "#include <g4main/PHG4Detector.h>\n";
    print F "\n";

    print F "#include <Geant4/G4Box.hh>\n";
    print F "#include <Geant4/G4Color.hh>\n";
    print F "#include <Geant4/G4LogicalVolume.hh>\n";
    print F "#include <Geant4/G4Material.hh>\n";
    print F "#include <Geant4/G4PVPlacement.hh>\n";
    print F "#include <Geant4/G4SystemOfUnits.hh>\n";
    print F "#include <Geant4/G4VisAttributes.hh>\n";
    print F "\n";

    print F "#include <cmath>\n";
    print F "#include <iostream>\n";
    print F "\n";

    print F "class G4VSolid;\n";
    print F "class PHCompositeNode;\n";
    print F "\n";

    print F "using namespace std;\n";
    print F "\n";

    print F "//____________________________________________________________________________..\n";
    print F "$detectorclassname\:\:$detectorclassname(PHG4Subsystem *subsys,\n";
    print F "                                         PHCompositeNode *Node,\n";
    print F "                                         PHParameters *parameters,\n";
    print F "                                         const std::string &dnam)\n";
    print F "  : PHG4Detector(subsys, Node, dnam)\n";
    print F "  , m_Params(parameters)\n";
    print F "{\n";
    print F "}\n";
    print F "\n";

    print F "//_______________________________________________________________\n";
    print F "int $detectorclassname\:\:IsInDetector(G4VPhysicalVolume *volume) const\n";
    print F "{\n";
    print F "  set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);\n";
    print F "  if (iter != m_PhysicalVolumesSet.end())\n";
    print F "  {\n";
    print F "    return 1;\n";
    print F "  }\n";
    print F "  return 0;\n";
    print F "}\n";
    print F "\n";

    print F "//_______________________________________________________________\n";
    print F "void $detectorclassname\:\:ConstructMe(G4LogicalVolume *logicWorld)\n";
    print F "{\n";
    print F " //begin implement your own here://\n";
    print F " // Do not forget to multiply the parameters with their respective CLHEP/G4 unit !\n";
    print F "  double xdim = m_Params->get_double_param(\"size_x\") * cm;\n";
    print F "  double ydim = m_Params->get_double_param(\"size_y\") * cm;\n";
    print F "  double zdim = m_Params->get_double_param(\"size_z\") * cm;\n";
    print F "  G4VSolid *solidbox = new G4Box(\"${detectorname}Solid\", xdim / 2., ydim / 2., zdim / 2.);\n";
    print F "  G4LogicalVolume *logical = new G4LogicalVolume(solidbox, G4Material::GetMaterial(m_Params->get_string_param(\"material\")), \"${detectorname}Logical\");\n";
    print F "\n";

    print F "  G4VisAttributes *vis = new G4VisAttributes(G4Color(G4Colour::Grey()));  // grey is good to see the tracks in the display\n";
    print F "  vis->SetForceSolid(true);\n";
    print F "  logical->SetVisAttributes(vis);\n";
    print F "  G4RotationMatrix *rotm = new G4RotationMatrix();\n";
    print F "  rotm->rotateX(m_Params->get_double_param(\"rot_x\") * deg);\n";
    print F "  rotm->rotateY(m_Params->get_double_param(\"rot_y\") * deg);\n";
    print F "  rotm->rotateZ(m_Params->get_double_param(\"rot_z\") * deg);\n";
    print F "\n";

    print F "  G4VPhysicalVolume *phy = new G4PVPlacement(\n";
    print F "      rotm,\n";
    print F "      G4ThreeVector(m_Params->get_double_param(\"place_x\") * cm,\n";
    print F "                    m_Params->get_double_param(\"place_y\") * cm,\n";
    print F "                    m_Params->get_double_param(\"place_z\") * cm),\n";
    print F "      logical, \"$detectorname\", logicWorld, 0, false, OverlapCheck());\n";
    print F "  // add it to the list of placed volumes so the IsInDetector method\n";
    print F "  // picks them up\n";
    print F "  m_PhysicalVolumesSet.insert(phy);\n";
    print F " //end implement your own here://\n";
    print F "  return;\n";
    print F "}\n";
    print F "\n";

    print F "//_______________________________________________________________\n";
    print F "void $detectorclassname\:\:Print(const std::string &what) const\n";
    print F "{\n";
    print F "  std::cout << \"$detectorname Detector:\" << std::endl;\n";
    print F "  if (what == \"ALL\" || what == \"VOLUME\")\n";
    print F "  {\n";
    print F "    std::cout << \"Version 0.1\" << std::endl;\n";
    print F "    std::cout << \"Parameters:\" << std::endl;\n";
    print F "    m_Params->Print();\n";
    print F "  }\n";
    print F "  return;\n";
    print F "}\n";

    close(F);
}

sub CreateDetectorInclude()
{
    my $file = shift;
    open(F,">$file");
    my $includeguard = uc(sprintf("%s_H",$detectorclassname));
    open(F,">$file");
    print F "// Tell emacs that this is a C++ source\n";
    print F "//  -*- C++ -*-.\n";
    print F "#ifndef $includeguard\n";
    print F "#define $includeguard\n";
    print F "\n";

    print F "#include <g4main/PHG4Detector.h>\n";
    print F "\n";

    print F "#include <set>\n";
    print F "#include <string>  // for string\n";
    print F "\n";

    print F "class G4LogicalVolume;\n";
    print F "class G4VPhysicalVolume;\n";
    print F "class PHCompositeNode;\n";
    print F "class PHG4Subsystem;\n";
    print F "class PHParameters;\n";
    print F "\n";

    print F "class $detectorclassname : public PHG4Detector\n";
    print F "{\n";
    print F " public:\n";
    print F "  //! constructor\n";
    print F "  $detectorclassname(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);\n";
    print F "\n";

    print F "  //! destructor\n";
    print F "  virtual ~$detectorclassname() {}\n";
    print F "\n";

    print F "  //! construct\n";
    print F "  void ConstructMe(G4LogicalVolume *world) override;\n";
    print F "\n";

    print F "  void Print(const std::string &what = \"ALL\") const override;\n";
    print F "\n";

    print F "  //!\@name volume accessors\n";
    print F "  //\@{\n";
    print F "  int IsInDetector(G4VPhysicalVolume *) const;\n";
    print F "  //\@}\n";
    print F "\n";

    print F "  void SuperDetector(const std::string &name) { m_SuperDetector = name; }\n";
    print F "  const std::string SuperDetector() const { return m_SuperDetector; }\n";
    print F "\n";

    print F " private:\n";
    print F "  PHParameters *m_Params;\n";
    print F "\n";

    print F "  // active volumes\n";
    print F "  std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;\n";
    print F "\n";

    print F "  std::string m_SuperDetector;\n";
    print F "};\n";
    print F "\n";

    print F "#endif // $includeguard\n";
    close(F);
}


sub CreateSubsystemInclude()
{
    my $includefile = shift;
    my $includeguard = uc(sprintf("%s_H",$subsysclassname));
    open(F,">$includefile");
    print F "// Tell emacs that this is a C++ source\n";
    print F "//  -*- C++ -*-.\n";
    print F "#ifndef $includeguard\n";
    print F "#define $includeguard\n";
    print F "\n";

    print F "#include <g4detectors/PHG4DetectorSubsystem.h>\n";
    print F "\n";

    print F "class PHCompositeNode;\n";
    print F "class PHG4Detector;\n";
    print F "class $detectorclassname;\n";
    print F "class PHG4SteppingAction;\n";
    print F "\n";
    print F "/**\n";
    print F "   * \\brief Detector Subsystem module\n";
    print F "   *\n";
    print F "   * The detector is constructed and registered via $detectorclassname\n";
    print F "   *\n";
    print F "   *\n";
    print F "   * \\see $detectorclassname\n";
    print F "   * \\see $subsysclassname\n";
    print F "   *\n";
    print F "   */\n";
    print F "class $subsysclassname : public PHG4DetectorSubsystem\n";
    print F "{\n";
    print F " public:\n";
    print F "  //! constructor\n";
    print F "  $subsysclassname(const std::string& name = \"$detectorname\");\n";
    print F "\n";

    print F "  //! destructor\n";
    print F "  virtual ~$subsysclassname() {}\n";
    print F "\n";

    print F "  /*!\n";
    print F "  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST\n";
    print F "  */\n";
    print F "  int InitRunSubsystem(PHCompositeNode*) override;\n";
    print F "\n";

    print F "  //! event processing\n";
    print F "  /*!\n";
    print F "  get all relevant nodes from top nodes (namely hit list)\n";
    print F "  and pass that to the stepping action\n";
    print F "  */\n";
    print F "  int process_event(PHCompositeNode*) override;\n";
    print F "\n";

    print F "  //! accessors (reimplemented)\n";
    print F "  PHG4Detector* GetDetector() const override;\n";
    print F "\n";

    print F "  PHG4SteppingAction* GetSteppingAction() const override { return m_SteppingAction; }\n";
    print F "  //! Print info (from SubsysReco)\n";
    print F "  void Print(const std::string& what = \"ALL\") const override;\n";
    print F "\n";

    print F " protected:\n";
    print F "  // \\brief Set default parameter values\n";
    print F "  void SetDefaultParameters() override;\n";
    print F "\n";

    print F " private:\n";
    print F "  //! detector construction\n";
    print F "  /*! derives from PHG4Detector */\n";
    print F "  $detectorclassname  *m_Detector;\n";
    print F "\n";

    print F "  //! particle tracking \"stepping\" action\n";
    print F "  /*! derives from PHG4SteppingActions */\n";
    print F "  PHG4SteppingAction *m_SteppingAction;\n";
    print F "};\n";
    print F "\n";

    print F "#endif // $includeguard\n";
    close(F);
}

sub CreateSubsystemImplementation()
{
    my $file = shift;
    open(F,">$file");
    print F "//____________________________________________________________________________..\n";
    print F "//\n";
    print F "// This is the interface to the framework. You only need to define the parameters\n";
    print F "// you use for your detector in the SetDefaultParameters() method here\n";
    print F "// The place to do this is marked by //implement your own here//\n";
    print F "// The parameters have no units, they need to be converted in the\n";
    print F "// $detectorclassname\:\:ConstructMe() method\n";
    print F "// but the convention is as mentioned cm and deg\n";
    print F "//____________________________________________________________________________..\n";
    print F "//\n";
    print F "#include \"$subsystem_includefile\"\n";
    print F "\n";
    print F "#include \"$detector_includefile\"\n";
    print F "#include \"$steppingaction_includefile\"\n";
    print F "\n";

    print F "#include <phparameter/PHParameters.h>\n";
    print F "\n";

    print F "#include <g4main/PHG4HitContainer.h>\n";
    print F "#include <g4main/PHG4SteppingAction.h> \n";
    print F "\n";

    print F "#include <phool/PHCompositeNode.h>\n";
    print F "#include <phool/PHIODataNode.h>\n";
    print F "#include <phool/PHNode.h> \n";
    print F "#include <phool/PHNodeIterator.h>\n";
    print F "#include <phool/PHObject.h>\n";
    print F "#include <phool/getClass.h>\n";
    print F "\n";

    print F "using namespace std;\n";
    print F "\n";

    print F "//_______________________________________________________________________\n";
    print F "$subsysclassname\:\:$subsysclassname(const std::string &name)\n";
    print F "  : PHG4DetectorSubsystem(name)\n";
    print F "  , m_Detector(nullptr)\n";
    print F "  , m_SteppingAction(nullptr)\n";
    print F "{\n";
    print F "  // call base class method which will set up parameter infrastructure\n";
    print F "  // and call our SetDefaultParameters() method\n";
    print F "  InitializeParameters();\n";
    print F "}\n";

    print F "//_______________________________________________________________________\n";
    print F "int $subsysclassname\:\:InitRunSubsystem(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  PHNodeIterator iter(topNode);\n";
    print F "  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(\"PHCompositeNode\", \"DST\"));\n";
    print F "  PHNodeIterator dstIter(dstNode);\n";
    print F "  if (GetParams()->get_int_param(\"active\"))\n";
    print F "  {\n";
    print F "    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst(\"PHCompositeNode\", Name()));\n";
    print F "    if (!DetNode)\n";
    print F "    {\n";
    print F "      DetNode = new PHCompositeNode(Name());\n";
    print F "      dstNode->addNode(DetNode);\n";
    print F "    }\n";
    print F "    std::string g4hitnodename = \"G4HIT_\" + Name();\n";
    print F "    PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(DetNode, g4hitnodename);\n";
    print F "    if (!g4_hits)\n";
    print F "    {\n";
    print F "      g4_hits = new PHG4HitContainer(g4hitnodename);\n";
    print F "      DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, g4hitnodename, \"PHObject\"));\n";
    print F "    }\n";
    print F "  }\n";
    print F "  // create detector\n";
    print F "  m_Detector = new $detectorclassname(this, topNode, GetParams(), Name());\n";
    print F "  m_Detector->OverlapCheck(CheckOverlap());\n";
    print F "  // create stepping action if detector is active\n";
    print F "  if (GetParams()->get_int_param(\"active\"))\n";
    print F "  {\n";
    print F "    m_SteppingAction = new $steppingclassname(m_Detector, GetParams());\n";
    print F "  }\n";
    print F "  return 0;\n";
    print F "}\n";

    print F "//_______________________________________________________________________\n";
    print F "int $subsysclassname\:\:process_event(PHCompositeNode *topNode)\n";
    print F "{\n";
    print F "  // pass top node to stepping action so that it gets\n";
    print F "  // relevant nodes needed internally\n";
    print F "  if (m_SteppingAction)\n";
    print F "  {\n";
    print F "    m_SteppingAction->SetInterfacePointers(topNode);\n";
    print F "  }\n";
    print F "  return 0;\n";
    print F "}\n";

    print F "//_______________________________________________________________________\n";
    print F "void $subsysclassname\:\:Print(const string &what) const\n";
    print F "{\n";
    print F "  if (m_Detector)\n";
    print F "  {\n";
    print F "    m_Detector->Print(what);\n";
    print F "  }\n";
    print F "  return;\n";
    print F "}\n";
    print F "\n";

    print F "//_______________________________________________________________________\n";
    print F "PHG4Detector *$subsysclassname\:\:GetDetector(void) const\n";
    print F "{\n";
    print F "  return m_Detector;\n";
    print F "}\n";
    print F "\n";

    print F "//_______________________________________________________________________\n";
    print F "void $subsysclassname\:\:SetDefaultParameters()\n";
    print F "{\n";
    print F "  // sizes are in cm\n";
    print F "  // angles are in deg\n";
    print F "  // units should be converted to G4 units when used\n";
    print F "  //implement your own here//\n";
    print F "  set_default_double_param(\"place_x\", 0.);\n";
    print F "  set_default_double_param(\"place_y\", 0.);\n";
    print F "  set_default_double_param(\"place_z\", 0.);\n";
    print F "  set_default_double_param(\"rot_x\", 0.);\n";
    print F "  set_default_double_param(\"rot_y\", 0.);\n";
    print F "  set_default_double_param(\"rot_z\", 0.);\n";
    print F "  set_default_double_param(\"size_x\", 20.);\n";
    print F "  set_default_double_param(\"size_y\", 20.);\n";
    print F "  set_default_double_param(\"size_z\", 20.);\n";
    print F "\n";

    print F "  set_default_string_param(\"material\", \"G4_Cu\");\n";
    print F "}\n";

    close(F);
}

sub CreateAutogen()
{
    open(F,">autogen.sh");
    print F "#!/bin/sh\n";
    print F "srcdir=`dirname \$0`\n";
    print F "test -z \"\$srcdir\" && srcdir=.\n";
    print F "\n";
    print F "(cd \$srcdir; aclocal -I \${OFFLINE_MAIN}/share;\\\n";
    print F "libtoolize --force; automake -a --add-missing; autoconf)\n";
    print F "\n";
    print F "\$srcdir/configure  \"\$\@\"\n";
    close(F);
    chmod 0755, "autogen.sh";
}

sub CreateConfigure()
{
    my $loname = lc($detectorname);
    open(F,">configure.ac");
    print F "AC_INIT($loname,[1.00])\n";
    print F "AC_CONFIG_SRCDIR([configure.ac])\n";
    print F "\n";

    print F "AM_INIT_AUTOMAKE\n";
    print F "AC_PROG_CXX(CC g++)\n";
    print F "\n";

    print F "LT_INIT([disable-static])\n";
    print F "\n";

    print F "dnl   no point in suppressing warnings people should \n";
    print F "dnl   at least see them, so here we go for g++: -Wall\n";
    print F "if test \$ac_cv_prog_gxx = yes; then\n";
    print F "   CXXFLAGS=\"\$CXXFLAGS -Wall -Werror\"\n";
    print F "fi\n";
    print F "\n";

    print F "AC_CONFIG_FILES([Makefile])\n";
    print F "AC_OUTPUT\n";
    close(F);
}

sub CreateMakefile()
{
    open(F,">Makefile.am");
    print F "AUTOMAKE_OPTIONS = foreign\n";
    print F "\n";

    print F "AM_CPPFLAGS = \\\n";
    print F "  -I\$(includedir) \\\n";
    print F "  -I\$(OFFLINE_MAIN)/include \\\n";
    print F "  -I\$(ROOTSYS)/include\\\n";
    print F "  -I\$(G4_MAIN)/include \n";
    print F "\n";

    print F "AM_LDFLAGS = \\\n";
    print F "  -L\$(libdir) \\\n";
    print F "  -L\$(OFFLINE_MAIN)/lib \\\n";
    print F "  -L\$(OFFLINE_MAIN)/lib64\n";
    print F "\n";

    print F "pkginclude_HEADERS = \\\n";
    print F "  $subsystem_includefile\n";
    print F "\n";

    print F "lib_LTLIBRARIES = \\\n";
    print F "  lib$detectorname.la\n";
    print F "\n";

    print F "lib${detectorname}_la_SOURCES = \\\n";
    print F "  $subsystem_implementfile\\\n";
    print F "  $detector_implementfile\\\n";
    print F "  $steppingaction_implementfile\n";
    print F "\n";

    print F "lib${detectorname}_la_LIBADD = \\\n";
    print F "  -lphool \\\n";
    print F "  -lSubsysReco\\\n";
    print F "  -lg4detectors\\\n";
    print F "  -lg4testbench \n";
    print F "\n";

    print F "BUILT_SOURCES = testexternals.cc\n";
    print F "\n";

    print F "noinst_PROGRAMS = \\\n";
    print F "  testexternals\n";
    print F "\n";

    print F "testexternals_SOURCES = testexternals.cc\n";
    print F "testexternals_LDADD   = lib$detectorname.la\n";
    print F "\n";

    print F "testexternals.cc:\n";
    print F "\techo \"//*** this is a generated file. Do not commit, do not edit\" > \$\@\n";
    print F "\techo \"int main()\" >> \$\@\n";
    print F "\techo \"{\" >> \$\@\n";
    print F "\techo \"  return 0;\" >> \$\@\n";
    print F "\techo \"}\" >> \$\@\n";
    print F "\n";

    print F "clean-local:\n";
    print F "\trm -f \$(BUILT_SOURCES)\n";
    close(F);
}

