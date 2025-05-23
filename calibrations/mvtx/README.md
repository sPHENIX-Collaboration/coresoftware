# MVTX Calibration

## Building the package

A package should be built before running the scripts, or remove any mention of these packages in the Fun4All macros so you can just process the data. 

There is an included script called `setup-condor.sh` will create condor files to submit jobs in `$pwd/condor`. This is an example of how to run the script:

```bash
export $MYINSTALL=<your install directory path>
source setup-condor.sh $MYINSTALL
```
Help can be found by running `setup-condor.sh` with the `-h` flag.
```bash
source setup-condor.sh -h
```

## What is in the package?

These modules are located in the `calibration/mvtx` directory. The modules are:

1. `mvtxcalib`
	- Contains classes and subsystem modules for doing a Fake hit rate analysis on raw `.prdf` files.

#### Required `coresoftware` dependencies:

These modules are located in the `offline/packages/mvtx` directory. The modules are:

- `MvtxPixelDefs`
	- Contains definitions and methods for creating pixel keys for each pixel in the MVTX detector. The pixel key is a 64-bit integer that uniquely identifies each pixel in the detector.
- `MvtxPixelMask`
	- Class which loads pixel masks from the conditions database and provides methods for masking pixels in hitmaps.
- `MvtxHitMap`
	- Class which creates hitmaps from raw data files. Hitmap is a vector of pairs of pixel keys and the number of hits in that pixel.


### `mvtxcalib` module

This module contains classes and subsystem modules for doing a Fake hit rate analysis on raw `.prdf` files. The module contains classes for loading pixel masks from the conditions database, creating hitmaps from raw data, and calculating the fake hit rate. 

- `MvtxFakeHitRate`
	- Main class for calculating the fake hit rate. Creates hitmaps from raw data files, masks pixels, and calculates the fake hit rate.
	- `SetOutputfile(std::string outputfile)`
		- Sets the output file for the fake hit rate analysis. The output file will contain the fake hit rate and the optimized pixel mask.
	- `SetMaxMaskedPixels(int max_masked_pixels)`
		- Sets the maximum number of pixels to mask in the optimization. The optimization will mask the `max_masked_pixels` pixels with the most hits in the hitmap.

- `MvtxHitMapWriter`
	- Writes the hitmap to a root file. The hitmap is a vector of pairs of pixel keys and the number of hits in that pixel. The hitmap is written to a root file with the name `mvtx_hit_map.root` by default. The root file contains a TTree with the hitmap data.
	- `SetOutputfile(std::string outputfile)`
		- Sets the output file for the hitmap analysis. The output file will contain the hitmap data.

## Running the Fake Hit Rate Analysis

The fake hit rate analysis is run using the `Fun4All_MvtxFHR.C` macro in the`<SPHENIXMACROS>/calibrations/mvtx` directory. This macro will run the fake hit rate analysis on raw MVTX data files and save the results to an output file. The macro takes the following arguments:
	
- `nEvents` : Number of events to process
- `run_number` : Run number of the data
- `trigger_rate_kHz` : Trigger rate in kHz
- `output_name` : Output file name
- `selected_layer` : Selected layer to analyze (optional, -1 for all layers)
- `trigger_guard_output_name` : Output file name for trigger guard analysis (optional) 

The macro can be run using the following command:

```bash
root -l -b -q 'Fun4All_MvtxFHR.C(10000, 42641, 44, "output.root", "-1", "trigger_guard_output.root")'
```

The macro will initialze a Fun4All server,  then load raw MVTX data files for the given run number using the `RawDataManager` class. The macro will then create a `MvtxFakeHitRate` object and set the output file name and the maximum number of pixels to mask. The macro will then run the fake hit rate analysis and save the results to the output file. The output file will contain the fake hit rate and the optimized pixel mask.

## Generating Pixel Masks

To generate a json file with the pixel mask, you can use the `MakeJsonHotPixelMap.C` macro in the `<SPHENIXMACROS>/calibrations/mvtx` directory. This macro will read the output file from the fake hit rate analysis and generate a json file with the pixel mask. The macro takes the following arguments:

- `calibration_file` : Output file from the fake hit rate analysis
- `target_threshold` : Target threshold for fake hit rate

The macro can be run using the following command:

```bash
root -l -b -q 'MakeJsonHotPixelMap.C("output.root", 10e-8)'
```

The macro will read the output file from the fake hit rate analysis and generate a json file with the pixel mask. The json file will contain the pixel mask and the fake hit rate for each pixel. The macro will also print the number of pixels masked and the fake hit rate threshold.


## Questions

This is a work in progress and likely has some bugs. If you have any questions or need help, please contact me at
` <tmengel@bnl.gov>`.
```



