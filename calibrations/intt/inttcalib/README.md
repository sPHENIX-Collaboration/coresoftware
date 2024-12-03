# InttCalib Module

The `InttCalib` class is a main of the calibration module for the INTT detector. It inherits from `SubsysReco` and provides methods for initializing, processing events, and finalizing the calibration process.

## Class: InttCalib

### Constructor
- `InttCalib(std::string const& name = "InttCalib")`
  - Initializes the `InttCalib` object with the given name.

### Destructor
- `~InttCalib()`
  - Cleans up dynamically allocated memory.

### Methods

#### Initialization and Run Control
- `int InitRun(PHCompositeNode* topNode) override`
  - Initializes the run, setting up necessary data structures and loading maps from the CDB.

- `int process_event(PHCompositeNode* topNode) :woverride`
  - Processes each event, updating hitmaps and handling BCO offsets.

- `int EndRun(int const run_number) override`
  - Finalizes the run, generating hot maps and BCO maps.

#### OutPutFile Setup
- `void SetHotMapCdbFile(std::string const& file)`
- `void SetHotMapPngFile(std::string const& file)`
- `void SetBcoMapCdbFile(std::string const& file)`
- `void SetBcoMapPngFile(std::string const& file)`

#### Configuration
- `void SetStreamingMode(bool mode)`
  - Set Streaming mode. Especially important for BCO Calibration in pp runs
- `void SetBcoMaximumEvent(int mext)`
  - Set Maxmium event to produce the BCO maps.
- `void SetRunNumber(int runnum)`
  - RunNumber Setting
- `void SetDoFeebyFee(bool in)`
  - true : Doing fitting fee by fee for accurate cold channel determination  (Default : false)
- `void SetColdSigmaCut(double in)`
  - Sigma cut for cold channel
- `void SetHotSigmaCut(double in)` 
  - Sigma cut for hot channel
- `void SetppMode(bool mode)`
  - Not used. Keeping for possibility of future implementation. 

#### Hot Map Generation / FELIX server by server
- `int ConfigureHotMap_v3()`
- `int MakeHotMapCdb_v3()`
- `int MakeHotMapPng_v3()`

#### Hot Map generation / Fee by fee (fee = 1 half ladder)
- `int ConfigureHotMap_fee()`
- `int MakeHotMapCdb_fee()`
- `int MakeHotMapROOT_fee()`

#### BCO Map Generation
- `int ConfigureBcoMap()`
- `int MakeBcoMapCdb()`
- `int MakeBcoMapPng()`

#### Histogram Configuration
- `int ConfigureHist_v3(TH1D*& hist, TF1*& fit, double maxbin, std::map<double, int> const& hitrate_map, std::string const& name, std::string const& title)`

#### Old versions (not used)
- `int ConfigureHotMap_v2()`
- `int MakeHotMapCdb_v2()`
- `int MakeHotMapPng_v2()`
- `int ConfigureHotMap()`
- `int MakeHotMapCdb()`
- `int MakeHotMapPng()`
- `int ConfigureHist(TH1D*& hist, TF1*& fit, std::map<double, int> const& hitrate_map, std::string const& name, std::string const& title)`
- `int ConfigureHist_v2(TH1D*& hist, TF1*& fit, std::map<double, int> const& hitrate_map, std::string const& name, std::string const& title)`



#### Utility Methods
- `int adjust_hitrate(InttMap::Offline_s const& ofl, double& hitrate) const`
- `int GetIndex(InttMap::RawData_s const& raw, InttMap::Offline_s const& ofl) const`
- `int GetFeeIndex(InttMap::RawData_s const& raw, InttMap::Offline_s const& ofl) const`
- `std::pair<double, double> CalculateStandardDeviation(const std::vector<int>& data)`
- `Color_t GetFeeColor(int fee) const`

#### Debugging and Data Management
- `void Debug()`
- `int SaveHitrates()`
- `int LoadHitrates()`

This class provides a comprehensive set of methods for calibrating the INTT detector, including handling hot and cold channels, generating maps, and configuring histograms.
