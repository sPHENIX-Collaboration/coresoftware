# TpcCalibrationQA

Compiled postprocessing of existing TpcSeedsQA and TPC timing histograms.

`TpcCalibrationQA` configures and runs histogram analysis. It also exposes
cache and collection entry points for the production driver. The internal
classes `TpcDistortionQA` and `TpcDriftQA` implement distortion and timing
analysis. Input validation and ROOT output/rendering live in separate sources.

`TpcConditionsReference` exports full-run CDB conditions and a PHGarfield
prediction. Its prediction uses the `TpcPHGarfieldReference` SubsysReco for
one-time run initialization. It is not an event histogram producer. Run each
export in a separate ROOT process from active reconstruction.

Build out of source with the usual sPHENIX environment:

```bash
/path/to/coresoftware/offline/QA/TpcCalibrationQA/autogen.sh --prefix="$MYINSTALL"
make -j4
make install
```

Installed public headers are under `tpccalibrationqa/`; the library is
`libtpccalibrationqa.so`. Use the matching macro and Python driver in
`macros/calibrations/tpc`.

The output schema and numerical algorithms are preserved from the original
ROOT macros. Physical-side and legacy PCA-z populations stay distinct.
Current/base-kEff snapshots do not certify historically applied corrections.

Histogram analysis has no event input manager. Concurrent jobs must use
separate processes; ROOT state and the output index are process-local.
