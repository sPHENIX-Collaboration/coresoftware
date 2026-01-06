

#include "AnnularFieldSim.h"

#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TFileCollection.h>
#include <TFileInfo.h>
#include <TH3.h>
#include <THashList.h>
#include <TTree.h>

#include <format>
#include <string>

R__LOAD_LIBRARY(libfieldsim.so)

std::string field_string;
std::string lookup_string;

AnnularFieldSim *SetupDefaultSphenixTpc(bool twinMe=false, bool useSpacecharge=true, float xshift=0, float yshift=0, float zshift=0);
AnnularFieldSim *SetupDigitalCurrentSphenixTpc(bool twinMe=false, bool useSpacecharge=true);
void TestSpotDistortion(AnnularFieldSim *t);
void SurveyFiles(TFileCollection* filelist);

  
void generate_distortion_map(const char *inputname, const char* gainName, const char *outputname, const char *ibfName, const char *primName, bool hasSpacecharge=true, bool isAdc=false, int  /*nSteps*/=500, bool scanSteps=false, float xshift=0, float yshift=0, float zshift=0){
  std::cout << "generating single distortion map.  Caution:  This is vastly less efficient than re-using the tpc model once it is set up" << std::endl;
 
  bool hasTwin=true; //this flag prompts the code to build both a positive-half and a negative-half for the TPC, reusing as much of the calculations as possible.  It is more efficient to 'twin' one half of the TPC than to recalculate/store the greens functions for both.
  TString gainHistName[2]={"hIbfGain_posz","hIbfGFain_negz"};

  //and some parameters of the files we're loading:
  bool usesChargeDensity=false; //true if source hists contain charge density per bin.  False if hists are charge per bin.
  float tpc_chargescale=1.6e-19;//Coulombs per bin unit.
  float spacecharge_cm_per_axis_unit=0.1;//cm per histogram axis unit (mm), matching the MDC2 sample from Evgeny.
  
  //file names we'll be filling as we go:
  TFile *infile;
  TFile *gainfile;
 
  TString sourcefilename=inputname;
  std::string outputfilename=outputname;

  //now build the time-consuming part:
  AnnularFieldSim *tpc;
    tpc=SetupDefaultSphenixTpc(hasTwin,hasSpacecharge, xshift, yshift, zshift);//loads the lookup, fields, etc.
 
  //and the location to plot the fieldslices about:
 TVector3 pos=0.5*(tpc->GetOuterEdge()+tpc->GetInnerEdge());;
  pos.SetPhi(3.14159);
  
  infile=TFile::Open(sourcefilename.Data(),"READ");


  //the total charge is prim + IBF
  //if we are doing ADCs, though, we only read the one.
  TH3* hCharge=(TH3*)(infile->Get(ibfName));
  if (!isAdc){
    hCharge->Add((TH3*)(infile->Get(primName)));
  }   
    //hCharge->Scale(70);//Scaleing the histogram spacecharge by 100 times

  std::string chargestring;
	       
  //load the spacecharge into the distortion map generator:
  //  void load_spacecharge(TH3F *hist, float zoffset, float chargescale, float cmscale, bool isChargeDensity);
  if (!isAdc){
    chargestring=std::format("{}:({}+{})",inputname,ibfName,primName);
    tpc->load_spacecharge(hCharge,0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity, chargestring);
    if (hasTwin) { tpc->twin->load_spacecharge(hCharge,0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity);
}
  }
  if (isAdc){ //load digital current using the scaling:
    gainfile=TFile::Open(gainName,"READ");
    TH2* hGain[2];
    hGain[0]=(TH2*)(gainfile->Get(gainHistName[0]));
    chargestring=std::format("{}:(dc:{} g:{}:{})",inputname,ibfName,gainName,gainHistName[0].Data());
    tpc->load_digital_current(hCharge,hGain[0],tpc_chargescale,spacecharge_cm_per_axis_unit,chargestring);
    if (hasTwin) {
      hGain[1]=(TH2*)(gainfile->Get(gainHistName[1]));
      tpc->twin->load_digital_current(hCharge,hGain[1],tpc_chargescale,spacecharge_cm_per_axis_unit,chargestring);
    }
  }
  //build the electric fieldmap from the chargemap
  tpc->populate_fieldmap();
  if (hasTwin) {  tpc->twin->populate_fieldmap();
}

  //build the distortion maps from the fieldmaps and save it to the output filename.
  if (scanSteps){
  for(int i=0;i<10;i++){
    std::string study_filestring=std::format("study_file_changinginterval.steps{}.hist.root",50*(i+1));
  tpc->GenerateSeparateDistortionMaps(study_filestring,50*(i+1),1,1,1,1,true);
  //tpc->GenerateSeparateDistortionMaps(outputfilename.Data(),1,1,1,1,false);
}
}
  else{
   tpc->GenerateSeparateDistortionMaps(outputfilename,450,1,1,1,1,false);
}

  std::cout << "distortions mapped." << std::endl;
  tpc->PlotFieldSlices(outputfilename,pos, 'E'); //plot the electric field
  tpc->PlotFieldSlices(outputfilename,pos,'B'); //plot the magnetic field
  std::cout << "fieldslices plotted." << std::endl;
  
  infile->Close();
  std::cout << "input closed.  All done.  Anything else is root's problem." << std::endl;

  return;
  
}

  
void generate_distortion_map(const char * inputpattern="./evgeny_apr/Smooth*.root", const char *outputfilebase="./apr07_maps/apr07", bool hasSpacecharge=true, bool isDigitalCurrent=false, int nSteps=500){
  

  int maxmaps=10;
//  int maxmapsperfile=2;

  
  bool hasTwin=true; //this flag prompts the code to build both a positive-half and a negative-half for the TPC, reusing as much of the calculations as possible.  It is more efficient to 'twin' one half of the TPC than to recalculate/store the greens functions for both.
  //bool hasSpacecharge=true;

  //and some parameters of the files we're loading:
  bool usesChargeDensity=false; //true if source hists contain charge density per bin.  False if hists are charge per bin.
  float tpc_chargescale=1.6e-19;//Coulombs per bin unit.
  float spacecharge_cm_per_axis_unit=0.1;//cm per histogram axis unit, matching the MDC2 sample from Evgeny.
  
  //file names we'll be filling as we go:
  TFile *infile;
 
  TString sourcefilename;
  std::string outputfilename;

  
  //find all files that match the input string (we don't need this yet, but doing it before building the fieldsim saves time if there's an empty directory or something.
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  filelist->Print();
  std::cout << "Using pattern \"" << inputpattern << "\", found " <<  filelist->GetList()->GetEntries()
	    << " files to read, eg : " << ((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetUrl() << std::endl;


  // SurveyFiles( filelist);

  //now build the time-consuming part:
  //AnnularFieldSim *tpc=SetupDefaultSphenixTpc(hasTwin,hasSpacecharge);//loads the lookup, fields, etc.
  AnnularFieldSim *tpc;


  //previously I flagged digital current and used a different rossegger lookup table.
  //now I just resample the histogram in that event.  There is a ~hard problem there to get right:
  //in each direction, if the source hist is binned more finely than the dest hist, I need to sum the total charge
  //if the source hist is binned more coarsely, I need to interpolate the charge density and multiply by the dest cell volume
  //but really, I need to differentiate between the source geometry and field point geometry.  task for another time...
  // if (isDigitalCurrent){
  //    tpc=SetupDigitalCurrentSphenixTpc(hasTwin,hasSpacecharge);//loads the lookup, fields, etc.
  //  }else {
    tpc=SetupDefaultSphenixTpc(hasTwin,hasSpacecharge);//loads the lookup, fields, etc.
    //  }
  //and the location to plot the fieldslices about:
 TVector3 pos=0.5*(tpc->GetOuterEdge()+tpc->GetInnerEdge());;
  pos.SetPhi(3.14159);


  


  double totalQ=0;
  
  for (int i=0;i<filelist->GetNFiles();i++){
   //for each file, find all histograms in that file.
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetFile();//gross
    infile=TFile::Open(sourcefilename.Data(),"READ");
    TList *keys=infile->GetListOfKeys();
    keys->Print();
    int nKeys=infile->GetNkeys();

    for (int j=0;j<nKeys;j++){
      TObject *tobj=infile->Get(keys->At(j)->GetName());
      //if this isn't a 3d histogram, skip it:
      bool isHist=tobj->InheritsFrom("TH3");
      if (!isHist) { continue;
}
      TString objname=tobj->GetName();
      if (objname.Contains("IBF")) { continue; //this is an IBF-only map we don't want.
}
      //assume this histogram is a charge map.
      //load just the averages:
      if (hasSpacecharge){
	if (isDigitalCurrent){
	  tpc->load_and_resample_spacecharge(8,36,40,sourcefilename.Data(),tobj->GetName(),0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity);
	  if (hasTwin) { tpc->twin->load_and_resample_spacecharge(8,36,40,sourcefilename.Data(),tobj->GetName(),0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity);
}

	}

	//load the spacecharge into the distortion map generator:
	tpc->load_spacecharge(sourcefilename.Data(),tobj->GetName(),0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity);
	if (hasTwin)
	{
	  tpc->twin->load_spacecharge(sourcefilename.Data(),tobj->GetName(),0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity);
	}

	//or load just the average:
	//tpc->load_spacecharge(sourcefilename.Data(),"h_Charge_0",0,tpc_chargescale,spacecharge_cm_per_axis_unit, usesChargeDensity);
	//printf("Sanity check:  Q has %d elements and dim=%d\n",tpc->q->Length(), tpc->q->dim);
	//totalQ=0;
	//for (int k=0;k<tpc->q->Length();k++){
	//  totalQ+=*(tpc->q->GetFlat(k));
	//}
	std::cout << std::format("Sanity check:  Total Q in reported region is {:E} C",totalQ) << std::endl;
      }
      tpc->populate_fieldmap();
      if (hasTwin)
      {
	tpc->twin->populate_fieldmap();
      }
      if (sourcefilename.Contains("verage")){ //this is an IBF map we don't want.
	outputfilename=std::format("{}.average.{}.{}",outputfilebase,field_string,lookup_string);
      } else {
	outputfilename=std::format("{}.file{}.{}.{}.{}",outputfilebase,i,tobj->GetName(),field_string,lookup_string);
      }
      std::cout << sourcefilename.Data() << " file has " << tobj->GetName() << " hist.  field="
		<< field_string << ", lookup=" << lookup_string << ". no scaling" << std::endl;

      //TestSpotDistortion(tpc);
 
      //tpc->GenerateDistortionMaps(outputfilename,2,2,2,1,true);


      tpc->GenerateSeparateDistortionMaps(outputfilename,2,2,2,1,nSteps,true);

      std::cout << "distortions mapped." << std::endl;
      tpc->PlotFieldSlices(outputfilename,pos);
      tpc->PlotFieldSlices(outputfilename,pos,'B');
      std::cout << "fieldslices plotted." << std::endl;     
      std::cout << "obj " << j << ": getname: " << tobj->GetName() << " inherits from TH3D:"
		<< tobj->InheritsFrom("TH3") << std::endl;

      //break; //rcc temp -- uncomment this to process one hist per file.
      if (i>maxmaps)
      {
	return;
      }
    }
      infile->Close();
  }

  return;
  
}

void TestSpotDistortion(AnnularFieldSim *t){
       TVector3 dummy(20.9034,-2.3553,-103.4712);
      float dummydest=-103.4752;
      t->twin->GetStepDistortion(dummydest,dummy,true,false);
      dummy.SetZ(dummy.Z()*-1);
      t->GetStepDistortion(-dummydest,dummy,true,false);
      return;
}

AnnularFieldSim *SetupDefaultSphenixTpc(bool twinMe, bool useSpacecharge, float xshift, float yshift, float zshift){
  //step1:  specify the sPHENIX space charge model parameters
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
//  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=105.5;
  const float tpc_cmVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  //const float tpc_magField=0.5;//T -- The old value used in carlos's studies.
  //const float tpc_driftVel=4.0*1e6;//cm per s  -- used in carlos's studies
  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  const float tpc_magField=1.4;//T -- 2019 nominal value
  const char detgeoname[]="sphenix";


  //step 2: specify the parameters of the field simulation.  Larger numbers of
  // bins will rapidly increase the memory footprint and compute times.  There
  // are some ways to mitigate this by setting a small region of interest, or a
  // more parsimonious lookup strategy, specified when AnnularFieldSim() is
  // actually constructed below.
  int nr=26;//10;//24;//159;//159 nominal
  int nr_roi_min=0;
  int nr_roi=nr;//10;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi=40;//38;//360;//360 nominal
  int nphi_roi_min=0;
  int nphi_roi=nphi;//38;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  int nz=40;//62;//62 nominal
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;

  bool realB=true;
  bool realE=true;
  

  
  //step 3:  create the fieldsim object.  different choices of the last few arguments will change how it
  // builds the lookup table spatially, and how it loads the spacecharge.  The various start-up steps
  // are exposed here so they can be timed in the macro.
  
  AnnularFieldSim *tpc;
  if (useSpacecharge){
    tpc=    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
				 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
  }else{
    tpc=    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
				 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::NoSpacecharge);
  }
    
  tpc->UpdateEveryN(10);//show reports every 10%.

    //load the field maps, either flat or actual maps
  tpc->setFlatFields(tpc_magField,tpc_cmVolt/tpc_z);
  field_string = std::format("flat_B{:.1f}_E{:.1f}",tpc_magField,tpc_cmVolt/tpc_z);

  if (realE){
    tpc->loadEfield("/sphenix/user/rcorliss/field/externalEfield.ttree.root","fTree");
    field_string = std::format("realE_B{:.1f}_E{:.1f}",tpc_magField,tpc_cmVolt/tpc_z);
  }
   if (realB){
    tpc->load3dBfield("/sphenix/user/rcorliss/field/sphenix3dmaprhophiz.root","fieldmap",1,-1.4/1.5, xshift, yshift,zshift);
        //tpc->loadBfield("sPHENIX.2d.root","fieldmap");
    field_string = std::format("realB_B{:.1f}_E{:.1f}_x{:.1f}_y{:.1f}_z{:.1f}",tpc_magField,tpc_cmVolt/tpc_z,xshift, yshift,zshift);
  } 
  if (realE && realB){
    field_string = std::format("real_B{:.1f}_E{:.1f}",tpc_magField,tpc_cmVolt/tpc_z);
  }
  std::cout << "set fields." << std::endl;




  //load the greens functions:
  lookup_string = std::format("ross_phi1_{}_phislice_lookup_r{}xp{}xz{}",detgeoname,nr,nphi,nz);
  //sprintf(lookupFilename,"%s.root",lookup_string);
  std::string lookupFilename = std::format("/sphenix/user/rcorliss/rossegger/{}.root",lookup_string); //hardcoded for racf
  TFile *fileptr=TFile::Open(lookupFilename.c_str(),"READ");

  if (!fileptr){ //generate the lookuptable
    //to use the full rossegger terms instead of trivial free-space greens functions, uncomment the line below:
    tpc->load_rossegger();
    std::cout << "loaded rossegger greens functions." << std::endl;
    tpc->populate_lookup();
    tpc->save_phislice_lookup(lookupFilename);
  } else{ //load it from a file
    fileptr->Close();
    tpc->load_phislice_lookup(lookupFilename);
  }

  std::cout << "populated lookup." << std::endl;






    
  //make our twin:
  if(twinMe){
    AnnularFieldSim *twin;
      if (useSpacecharge){
	twin=      new  AnnularFieldSim(tpc_rmin,tpc_rmax,-tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
      } else{
	twin=      new  AnnularFieldSim(tpc_rmin,tpc_rmax,-tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::NoSpacecharge);
      }    twin->UpdateEveryN(10);//show reports every 10%.

    //same magnetic field, opposite electric field
    twin->setFlatFields(tpc_magField,-tpc_cmVolt/tpc_z);
    //sprintf(field_string,"flat_B%2.1f_E%2.1f",tpc_magField,tpc_cmVolt/tpc_z);
    //twin->loadBfield("sPHENIX.2d.root","fieldmap");
    twin->load3dBfield("/sphenix/user/rcorliss/rossegger/sphenix3dmaprhophiz.root","fieldmap",1,-1.4/1.5, xshift, yshift, zshift);
    twin->loadEfield("/sphenix/user/rcorliss/rossegger/externalEfield.ttree.root","fTree",-1);//final '-1' tells it to flip z and the field z coordinate. r coordinate doesn't change, and we assume phi won't either, though the latter is less true.




    //borrow the greens functions:
    twin->borrow_rossegger(tpc->green,tpc_z);//use the original's green's functions, shift our internal coordinates by tpc_z when querying those functions.
    twin->borrow_epartial_from(tpc,tpc_z);//use the original's epartial.  Note that those values ought to be symmetric about z, and since our boundary conditions are translated along with our coordinates, they're completely unchanged.

    tpc->set_twin(twin);
  }

  return tpc;
}
  

AnnularFieldSim *SetupDigitalCurrentSphenixTpc(bool twinMe, bool useSpacecharge){
  //step1:  specify the sPHENIX space charge model parameters
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
//  float tpc_deltar=tpc_rmax-tpc_rmin;
  const float tpc_z=105.5;
  const float tpc_cmVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  //const float tpc_magField=0.5;//T -- The old value used in carlos's studies.
  //const float tpc_driftVel=4.0*1e6;//cm per s  -- used in carlos's studies
  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  const float tpc_magField=1.4;//T -- 2019 nominal value
  const char detgeoname[]="sphenix";
  
  //step 2: specify the parameters of the field simulation.  Larger numbers of
  // bins will rapidly increase the memory footprint and compute times.  There
  // are some ways to mitigate this by setting a small region of interest, or a
  // more parsimonious lookup strategy, specified when AnnularFieldSim() is
  // actually constructed below.
  
  int nr=8;//10;//24;//159;//159 nominal
  int nr_roi_min=0;
  int nr_roi=nr;//10;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi=3*12;//38;//360;//360 nominal
  int nphi_roi_min=0;
  int nphi_roi=nphi;//38;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  int nz=40;//62;//62 nominal
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;

  bool realB=true;
  bool realE=true;
  

  //step 3:  create the fieldsim object.  different choices of the last few arguments will change how it
  // builds the lookup table spatially, and how it loads the spacecharge.  The various start-up steps
  // are exposed here so they can be timed in the macro.
  AnnularFieldSim *tpc;
  if (useSpacecharge){
    tpc=    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
				 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
  }else{
    tpc=    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
				 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::NoSpacecharge);
  }
    
  tpc->UpdateEveryN(10);//show reports every 10%.

    //load the field maps, either flat or actual maps
  tpc->setFlatFields(tpc_magField,tpc_cmVolt/tpc_z);
  field_string = std::format("flat_B{:.1f}_E{:.1f}",tpc_magField,tpc_cmVolt/tpc_z);

  if (realE){
    tpc->loadEfield("externalEfield.ttree.root","fTree");
    field_string = std::format("realE_B{:.1f}_E{:.1f}",tpc_magField,tpc_cmVolt/tpc_z);
  }
   if (realB){
    tpc->load3dBfield("sphenix3dmaprhophiz.root","fieldmap",1,-1.4/1.5);
        //tpc->loadBfield("sPHENIX.2d.root","fieldmap");
    field_string = std::format("realB_B{:.1f}_E{:.1f}",tpc_magField,tpc_cmVolt/tpc_z);
  } 
  if (realE && realB){
    field_string = std::format("real_B{:.1f}_E{:.1f}",tpc_magField,tpc_cmVolt/tpc_z);
  }
  std::cout << "set fields." << std::endl;




  //load the greens functions:
  lookup_string = std::format("ross_phi1_{}_phislice_lookup_r{}xp{}xz{}",detgeoname,nr,nphi,nz);
  std::string lookupFilename = lookup_string + ".root";
  TFile *fileptr=TFile::Open(lookupFilename.c_str(),"READ");

  if (!fileptr){ //generate the lookuptable
    //to use the full rossegger terms instead of trivial free-space greens functions, uncomment the line below:
    tpc->load_rossegger();
    std::cout << "loaded rossegger greens functions." << std::endl;
    tpc->populate_lookup();
    tpc->save_phislice_lookup(lookupFilename);
  } else{ //load it from a file
    fileptr->Close();
    tpc->load_phislice_lookup(lookupFilename);
  }

  std::cout << "populated lookup." << std::endl;






    
  //make our twin:
  if(twinMe){
    AnnularFieldSim *twin;
      if (useSpacecharge){
	twin=      new  AnnularFieldSim(tpc_rmin,tpc_rmax,-tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
      } else{
	twin=      new  AnnularFieldSim(tpc_rmin,tpc_rmax,-tpc_z,
			   nr, nr_roi_min,nr_roi_max,1,2,
			   nphi,nphi_roi_min, nphi_roi_max,1,2,
			   nz, nz_roi_min, nz_roi_max,1,2,
			   tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::NoSpacecharge);
      }    twin->UpdateEveryN(10);//show reports every 10%.

    //same magnetic field, opposite electric field
    twin->setFlatFields(tpc_magField,-tpc_cmVolt/tpc_z);
    //sprintf(field_string,"flat_B%2.1f_E%2.1f",tpc_magField,tpc_cmVolt/tpc_z);
    //twin->loadBfield("sPHENIX.2d.root","fieldmap");
    twin->load3dBfield("sphenix3dmaprhophiz.root","fieldmap",1,-1.4/1.5);
    twin->loadEfield("externalEfield.ttree.root","fTree",-1);//final '-1' tells it to flip z and the field z coordinate. r coordinate doesn't change, and we assume phi won't either, though the latter is less true.




    //borrow the greens functions:
    twin->borrow_rossegger(tpc->green,tpc_z);//use the original's green's functions, shift our internal coordinates by tpc_z when querying those functions.
    twin->borrow_epartial_from(tpc,tpc_z);//use the original's epartial.  Note that those values ought to be symmetric about z, and since our boundary conditions are translated along with our coordinates, they're completely unchanged.

    tpc->set_twin(twin);
  }

  return tpc;
}
  
void  SurveyFiles(TFileCollection *filelist){
   TFile *infile;
 
  TString sourcefilename;
  TString outputfilename;
 //run a check of the files we'll be looking at:
  float integral_sum=0;
  int nhists=0;
  for (int i=0;i<filelist->GetNFiles();i++){
   //for each file, find all histograms in that file.
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetFile();//gross
    std::cout << "file " << i << ": " << sourcefilename << std::endl;
    infile=TFile::Open(sourcefilename.Data(),"READ");
    TList *keys=infile->GetListOfKeys();
    //keys->Print();
    int nKeys=infile->GetNkeys();

    for (int j=0;j<nKeys;j++){
      TObject *tobj=infile->Get(keys->At(j)->GetName());
      //if this isn't a 3d histogram, skip it:
      bool isHist=tobj->InheritsFrom("TH3");
      if (!isHist) { continue;
}
      TString objname=tobj->GetName();
      std::cout << " hist " << objname << " ";
      if (objname.Contains("IBF")) {
	std::cout << " is IBF only." << std::endl;
	continue; //this is an IBF map we don't want.
      }
      float integral=((TH3D*)tobj)->Integral();
      integral_sum+=integral;
      nhists+=1;
      std::cout << std::format(" will be used.  Total Q={:.3E} (ave={:.3E})",integral,integral_sum/nhists) << std::endl;
      //assume this histogram is a charge map.
      //load just the averages:
 
      //break; //rcc temp -- uncomment this to process one hist per file.
      //if (i>maxmaps) return;
    }
      infile->Close();
  }
  return;
}
