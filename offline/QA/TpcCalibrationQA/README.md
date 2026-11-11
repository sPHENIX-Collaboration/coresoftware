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

The output stores distortion profiles, timing edges and reference comparisons. Physical-side and legacy PCA-z populations stay distinct.
Current/base-kEff snapshots do not certify historically applied corrections.

Histogram analysis has no event input manager. Concurrent jobs must use
separate processes; ROOT state is process-local; each sample owns its output index.

## Drift input availability

Official upstream QA provides the `TpcRawHitQA` timing fallback. The optional
`h_TpcDriftVelocityQA_*` and `h_PHGarfieldDriftVelocityQA_*` inputs require
additional histogram-producing modules that are not part of upstream
coresoftware. This package reads those histograms when supplied; it does not
produce them. Missing ADC-weighted timing/z diagnostics are reported as
unavailable rather than inferred from thresholded sample counts.

Native calibration sums from the companion producer use `TH1D` (double
precision). Geometry and clock compatibility retain the existing absolute
`1e-6` tolerance in the respective cm/ns units. Float calibration-sum inputs
are not an established supported producer format.

## Batch recovery and central installation

New batch jobs use persistent `.lock` files with POSIX advisory locks. A killed
process releases its lock automatically; child ROOT processes retain it until
they exit. Do not delete `.lock` files: their presence does not mean a job is
running. An active job reports `Already running`. Legacy `.running` directories
are preserved: verify that the old job has stopped before removing one.

`run --allow-partial` permits per-run export failures and records/excludes
those runs during processing. Collection still requires a completed reference
run for every selected production group. Without the flag, failed exports
stop processing. Plotting flags control PDF creation only, not analysis ROOT
objects.

Central release installation additionally requires a registration in the
separate `sPHENIX-Collaboration/utilities` repository at
`utils/rebuild/packages.txt`, after the TPC QA entry. The source bundle includes
`utilities_packages.patch` for this companion change.
