# PHJewel — JEWEL 2.6.0 inside Fun4All

`PHJewel` is a SubsysReco that runs JEWEL in-process, the way `PHPythia8` wraps PYTHIA 8 and
`PHHerwig7` (ProjectHERWIG/Source/phherwig7) wraps Herwig 7. There is one Fun4All process per job,
with no FIFO and no HepMC on disk. The seed comes from `PHRandomSeed()`. Generated count, weight sum
and cross section go to `RUN/PHGenIntegral`. Truth triggers run inside the generation loop, so only
triggered events reach the node tree.

```
user card ──JewelCard (strict)──> canonical card ──> libjewel_<medium>.so (dlopen, RTLD_LOCAL)
                                                        │ init / initpythia / genevent  (JEWEL's own code)
                                                        │ HepMC2 lines in memory (JEWEL's own writer logic)
                       PHJewel <──HepMC::IO_GenEvent────┘
                          │ b [fm], PdfInfo, process id added
                          │ HepMCTruthTriggerCore (regenerate until pass)
                          └──> PHHepMCGenEventMap, RUN/PHGenIntegral ──> DST
```

| PHPythia8 | PHJewel |
|---|---|
| `Pythia::init()` | `init()` of JEWEL on the canonical card, seed into `MRPY(1)` |
| `Pythia::next()` | `genevent()` until JEWEL writes an event; isospin channel by `initpythia()` |
| `Pythia8ToHepMC` | JEWEL's `CONVERTTOHEPMC`, the same lines it writes to `HEPMCFILE`, read back by `HepMC::IO_GenEvent` |
| `info.sigmaGen()` | PYTHIA 6 `PARI(1)`, averaged over PYTHIA initialisations |
| `PHPy8JetTrigger` | `HepMCTruthTriggerCore`, the same code as PHHerwig7 and `hwfilter` |

## Validation status (2026-10-23)

| check | result |
|---|---|
| standalone binary rebuilt from the generated library sources vs the production `jewel-2.6.0-simple` (laptop), Au+Au 0-10% with recoils, 100 events | **byte-identical** HepMC |
| in-process API, Fortran writer on, same card, seed 7000 = NJOB 7 | **byte-identical** to production |
| in-process API, events from the memory buffer after the HepMC reader, Au+Au with recoils (100 events, 88k particles) and p+p vac (300 events) | **identical** within JEWEL's 9 printed digits (`tools/compare_hepmc.py`) |
| RCF, Fun4All + PHJewel → DST → read back, 30 Au+Au events | identical to the standalone binary from the same toolchain, and to the laptop production binary |
| RCF, jet trigger (R=0.4, 15 GeV, \|η\|<0.7) on p+p | 227 generated, 50 written, PHGenIntegral N/Σw/σ consistent |
| RCF, ISOCHANNEL XX with channel switching | PYTHIA re-initialised only on change; DST and PHGenIntegral fine |
| card rejection | the production cards `auau_0_10_recoil.dat` and `pp_recoil.dat` are refused with the reasons listed; a card with one of every mistake reports all of them in one pass |
| core libraries export only `jw_*` (RCF `nm -D`) | yes; PYTHIA 6 symbols and common blocks are private |

Not yet checked: Geant4 transport of the DST (PHG4 reading `PHHepMCGenEventMap`), embedding into a
background event, and long runs (memory, a JEWEL discard storm). These are the first things for Blair.

The RCF production binary `ProjectJEWEL/Code/jewel-2.6.0/jewel-2.6.0-simple` gives a **different event
stream** from the laptop binary for the same NJOB. Its PYTHIA, medium and subtraction objects date from
09-04, built by another compiler. JEWEL's showers are chaotic, so a floating-point difference changes
every event after the first. The samples are statistically equivalent, but they are not
event-for-event comparable. PHJewel on RCF matches the laptop build.

## Build

### RCF (sPHENIX)

```bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export MYINSTALL=/sphenix/tg/tg01/jets/nagle/install   # PHJewel and its trigger dependency live here
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
# HepMCTruthTriggerCore.h comes from ProjectHERWIG/Source/fun4all; it is installed into the same
# prefix (built in fun4all/build_tg), install it first if the prefix is new
cd phjewel && mkdir -p build && cd build && ../autogen.sh --prefix=$MYINSTALL && make -j4 install
```

This installs `libPHJewel.so` (the SubsysReco) and three cores, `libjewel_simple.so`, `libjewel_brick.so`
and `libjewel_vac.so`, each with its own private JEWEL plus PYTHIA 6.4.25 (`jewel/exports.map`). PHJewel
`dlopen`s the core named by `set_medium()` from `LD_LIBRARY_PATH`.

### Laptop (no Fun4All)

```bash
make -f Makefile.local       # lib/libjewel_*.dylib, bin/jewel-std-*, bin/jwgen_test
```

`bin/jewel-std-<medium>` is the standalone JEWEL rebuilt from the library sources, used as the
regression reference. `bin/jwgen_test` starts JEWEL exactly as PHJewel does and writes the events
after the HepMC reader:

```bash
PHJEWEL_CORE=lib/libjewel_simple.dylib bin/jwgen_test cards/auau200_0010_recoil.dat simple 100 7000 /tmp/x
```

### Library sources

`jewel/jewel-lib.f`, `jewel-main.f` and `medium-*.f` are **generated** from `../jewel-2.6.0` (the tree
the production binaries come from, with every patch in `../patches`):

```bash
make -f Makefile.local sources     # python3 tools/make_api_sources.py ../jewel-2.6.0 jewel
```

Never edit them by hand. The generator anchors every edit on exact text and stops if JEWEL changes
under it. The API hooks are inert unless `jw_init` set the API mode, so the standalone program built
from the same file is unchanged (see the first validation row). The hooks are:
- the parameter file is read from a common block instead of `getarg`;
- PHJewel owns the output files, NJOB and the seed;
- `out.ran` is not written;
- each `WRITE` in `CONVERTTOHEPMC` goes through a helper that writes the same line and also passes the numbers to the event buffer;
- file-name buffers are 512 characters instead of 80.

## Running

```bash
root -l -b -q 'Fun4All_PHJewel.C("cards/auau200_0010_recoil.dat","simple","out.root",1000,"jet",0.4,10.,0.7)'
```

```cpp
PHJewel *jw = new PHJewel();
jw->set_card("cards/auau200_0010_recoil.dat");   // main card
jw->set_medium("simple");                         // simple | brick | vac
jw->set_seed(0);              // 0: PHRandomSeed();  k*1000 reproduces a standalone NJOB k run
jw->set_output_stem("jw");    // jw.log, jw.card.dat, jw.medium.dat  (default PHJewel_<seed>)
jw->require_jet(); jw->set_jet_R(0.4); jw->set_jet_ptmin(10); jw->set_jet_etamax(0.7);
jw->set_max_trigger_attempts(100000);
// validation / extras
jw->set_fortran_hepmc_file("jw_fortran.hepmc");   // JEWEL's own writer, every JEWEL event
jw->set_evinfo_file("jw.evinfo");                 // EVINFOFILE incl. VIRTLOG columns
jw->set_isospin_block(1000);                      // ISOCHANNEL XX, see below
```

`nevents` in `se->run(n)` counts **events written** (passing the trigger), as for PHPythia8.

Only one PHJewel is allowed per process, because JEWEL keeps its state in Fortran common blocks. Do not
add PHPythia6 to the same process: the PYTHIA 6 symbols are isolated, but the Fortran I/O unit numbers
are shared by the runtime.

---

## Input cards

PHJewel reads **JEWEL's own card format**, so one file drives both PHJewel and standalone JEWEL. It reads
it **strictly**, because JEWEL's reader never fails. An unknown key, an unreadable value, a missing
medium file or a forgotten key all leave a built-in default in place, and those defaults are LHC
Pb+Pb (SQRTS 2760, MASS 208, TI 0.59, SIGMANN 6.2). A RHIC card with one typo runs to completion as a
different sample, and nothing in the output says so.

A run uses a **pair** of cards: the main card, plus the medium card named by `MEDIUMPARAMS` (not needed
for `vac`). The medium type is not in either card. `set_medium()` gives it, and the medium card must
contain exactly that medium's keys. A brick card given to `simple` fails on `TBRICK`.

### Grammar (the subset standalone JEWEL reads the same way)

| rule | why |
|---|---|
| `KEY value`: key in column 1, separated by **one space** | JEWEL splits at the first space; a key is everything before it |
| `#` comments **in column 1 only** | an indented `#` gives JEWEL an empty key and the line is ignored |
| no tabs, no DOS line endings | a tab ends up inside the key; `\r` ends up inside the value |
| key lines ≤ 80 characters (main), ≤ 100 (medium) | standalone JEWEL silently truncates longer lines; use relative paths |
| every key at most once | JEWEL takes the last one |
| integers without a decimal point (`NF 3`, not `3.`) | Fortran list-directed integer read fails on `3.` and keeps the default |
| reals as Fortran literals (`200.`, `1.5`, `1.d-3`) | `1,5` or `20GeV` fail the same way |
| logicals `T` / `F` | anything else fails the same way |
| paths: one token, nothing after it | JEWEL reads a path as the whole rest of the line, comment included |
| a trailing `# comment` is allowed after numbers and logicals | list-directed reads stop at the blank |
| ≤ 1000 lines | JEWEL reads the first 1000 |

**Relative paths are relative to the card file's directory.** Standalone JEWEL resolves them relative
to the working directory, so to run the same card standalone, run from the card's directory.

### Main card keys

**Required.** These have no default in PHJewel, because JEWEL's defaults are for the LHC:

| key | type | meaning |
|---|---|---|
| `PROCESS` | PPJJ PPYJ PPYQ PPYG PPZJ PPZQ PPZG PPWJ PPWQ PPWG PPDY | hard process (EEJJ refused) |
| `SQRTS` | real | √s_NN [GeV] |
| `PTMIN` | real ≥ 3 | p̂_T minimum [GeV]; JEWEL would raise anything lower to 3 without saying so |
| `PTMAX` | real | p̂_T maximum [GeV], `-1.` = none |
| `ETAMAX` | real > 0 | max \|η\| of ME partons and of the medium |
| `PDFSET` | int | LHAPDF ID of the proton PDF (13100 = CT14nlo) |
| `MASS`, `NPROTON` | int | nucleus A and Z: **1 1 for p+p** (JEWEL's default is 208 82) |
| `ISOCHANNEL` | XX PP PN NP NN | XX mixes nucleon pairs by Z/A; otherwise one pair only |
| `WEIGHTED` | T/F | p̂_T-weighted generation; with T, `WEXPO` is also required |
| `KEEPRECOILS` | T/F | keep medium recoils; with T, `WRITESCATCEN T` is required |
| `SPLITINTFILE`, `PDFFILE`, `XSECFILE` | path | JEWEL's tables; must exist (see Tables) |
| `MEDIUMPARAMS` | path | medium card (required for simple and brick, ignored for vac) |

**Optional.** JEWEL's built-in value applies if the key is absent:

| key | default | | key | default |
|---|---|---|---|---|
| `WEXPO` | 5. | | `HADROTYPE` | 0 (0 1) |
| `NF` | 3 (must equal medium NF) | | `PI0DEC` | T |
| `LAMBDAQCD` | 0.4 | | `ISRSCAT` | T |
| `LAMBDAQCDPS` | 0.4 | | `SHORTHEPMC` | T |
| `Q0` | 1.5 | | `COMPRESS` | T |
| `CHANNEL` | MUON (ELEC) | | `WRITESCATCEN` | F |
| `PDFALPHAS` | 0.118 | | `RECSOFTCUT` | 0. |
| `XSECSCALEFAC` | 1. | | `RECHARDCUT` | 5. |
| `MPI` | F | | `KINMODE` | 1 (0-3) |
| `MPIFSR` | T | | `RECMODE` | 0 (0-3; 2 turns recoils off) |
| `ANGORD` | 1 (0-5) | | `COHLENGTHFAC` | 1. |
| `AOFIRSTSPLIT` | F | | `PCOHREJ` | 1. |
| `SCATRECOIL` | F | | `COHSCAT` | T |
| `HADRO` | T (must be T) | | `RSUB`, `PUTONSHELL` | 1., 0 (no effect without DOSUBTRACTION) |

**Owned by PHJewel.** These are accepted but ignored, with a note in the log:

| key | replaced by |
|---|---|
| `NEVENT` | `se->run(n)` |
| `NJOB` | `set_seed()` / `PHRandomSeed()`; `MRPY(1) = seed`, so NJOB k ≡ seed 1000k |
| `LOGFILE` | `<stem>.log` |
| `HEPMCFILE` | the node tree (`set_fortran_hepmc_file` for validation) |
| `EVINFOFILE` | `set_evinfo_file` |
| `READRAN` | never |

**Refused values.** These are legal for JEWEL but wrong in front of a detector simulation:

| setting | why |
|---|---|
| `WRITEDUMMIES T` | adds one **status-1** π⁰ placeholder per thermal parton (E = 10⁻⁶ GeV along its direction, a ghost for constituent subtraction). There are hundreds per event; they carry no visible energy, but Geant and every truth-particle count (multiplicities, constituents, fragmentation) treat them as real particles. Recoil subtraction needs only `WRITESCATCEN T`. **The production cards `auau_*_recoil.dat` and `pp_recoil.dat` set it.** |
| `DOSUBTRACTION T` | rewrites the final state inside JEWEL; a DST needs the unsubtracted event |
| `KEEPRECOILS T` with `WRITESCATCEN F` | the recoil contribution could never be subtracted |
| `HADRO F` | partons in the final state |
| `PROCESS EEJJ` | not a hadron collider |
| `RECMODE 2` with `KEEPRECOILS T` | JEWEL silently turns recoils off |

### Medium cards

`simple` (Bjorken-expanding Glauber medium). Units: TAUI fm/c, temperatures GeV, **SIGMANN fm²
(42 mb = 4.2)**, centrality %.

| key | required | default | check |
|---|---|---|---|
| `TAUI`, `TI`, `TC` | yes | | TI > TC > 0; TI > 2 refused as "MeV" |
| `CENTRMIN`, `CENTRMAX` | yes | | 0 ≤ min < max ≤ 100 |
| `SIGMANN` | yes | | > 15 refused as "mb" |
| `NF` | yes | | = main-card NF |
| `MDFACTOR`, `MDSCALEFAC` | yes | | > 0 |
| `WOODSSAXON` | no | T | |
| `N0`, `D` | no | 0.17, 0.54 | > 0 |

Also refused: MASS < 2 with the simple medium.

`brick` (static uniform brick): `TBRICK` (GeV), `LBRICK` (fm), `NF`, `MDFACTOR`, `MDSCALEFAC`, all required.

`vac`: no medium card.

### What JEWEL is actually given

PHJewel writes `<stem>.card.dat` and `<stem>.medium.dat`, then hands those to JEWEL:
- keys that PHJewel owns are dropped, and paths are made absolute;
- every key you did not write appears as a `#default KEY value` comment, so the file is the complete record of the run;
- defaults are **not** written as active keys. JEWEL's built-in defaults are single-precision literals (`lqcd = 0.4` is 0.4000000059), so `LAMBDAQCD 0.4` would change the run in the 8th digit and break event identity with standalone JEWEL.

A 64-bit hash of the canonical pair goes into the log and into `PHGenIntegral::get_Description()`,
together with the medium, card path, seed and embedding id. After `init` PHJewel reads JEWEL's own
values back (A, Z, NF, √s, p̂_T range, η range, Q0, Λ, the recoil, subtraction and hadronisation
switches, process and isospin). It stops if any differ from its own reading of the card.

### Tables

`SPLITINTFILE`, `PDFFILE` and `XSECFILE` must exist. If one is missing, JEWEL integrates it, which is
slow, and writes it to the path, which is a race when parallel jobs share a directory.
`allow_table_generation()` permits that for a dedicated one-off run.

JEWEL checks only that a table exists, not what it was built with. The splitting integrals depend on
Q0 and the energy reach. The PDF and cross-section tables also depend on the medium's Debye-mass
range. Lookups outside a table extrapolate, and JEWEL counts them. **PHJewel prints the fractions at
End.** With the tables the production runs use (`pdfs_ct14nlo.dat`, `xsecs_ct14nlo.dat`), the RHIC
0-10 % simple medium extrapolates in about **11-12 % of PDF and cross-section lookups**. The brick and
vac media extrapolate in none. This is a property of the existing production samples too, and worth a
look independent of PHJewel.

### Cards shipped in `cards/`

| card | medium | |
|---|---|---|
| `pp200_vac.dat` | vac | p+p 200 GeV, p̂_T > 10, ISOCHANNEL PP |
| `auau200_vac.dat` | vac | Au+Au nucleon content (XX), no medium; JEWEL-only R_AA reference |
| `auau200_0010_recoil.dat` | simple | Au+Au 0-10 %, recoils on (WRITESCATCEN T) |
| `auau200_0010_norecoil.dat` | simple | the same, recoils off |
| `medium_auau200_0010.dat` | — | TAUI 0.8, TI 0.35, TC 0.17, SIGMANN 4.2, MDSCALEFAC 0.9 (arXiv:1212.1599 RHIC tune, as in the production studies) |
| `brick200_T30_L3.dat` + `brick_T30_L3.dat` | brick | T = 0.30 GeV, L = 3 fm |

`cards/tables/` holds symlinks to the production tables (laptop: `../jewel-2.6.0`, RCF:
`/sphenix/tg/tg01/jets/nagle/ProjectJEWEL/Code/jewel-2.6.0`).

---

## The event record

| | |
|---|---|
| status 1 | final state (what Geant transports and the trigger sees) |
| status 3 | thermal momenta of the scattering centres (`WRITESCATCEN T`): subtract them from jets at analysis level |
| status 4 | beam particles; also non-final partons |
| event number | JEWEL's NGOOD, which matches the EVINFOFILE line and the standalone HepMC file; after a trigger the numbers have gaps |
| weight | `weights()[0]` = PYTHIA `PARI(10)` |
| `GenCrossSection` | `PARI(1)` [pb] at that event, a running estimate for the current isospin channel |
| `HeavyIon::impact_parameter` | **b in fm** (simple), LBRICK (brick), -1 (vac); event-plane angle 0 because JEWEL lays b along +x. JEWEL itself writes the centrality in this slot, and PHJewel replaces it |
| `PdfInfo` | id1, id2 (`MSTI(15,16)`), x1, x2 (`PARI(33,34)`), Q = √`PARI(22)`; JEWEL leaves it empty |
| `signal_process_id` | PYTHIA `MSTI(1)` |
| vertices | all at the origin, as JEWEL writes them; the collision vertex comes from PHHepMCGenHelper's settings |

### Normalisation

`RUN/PHGenIntegral` holds `N_Generator_Accepted_Event` (events JEWEL wrote, before the trigger),
`N_Processed_Event` (events in the DST) and `Sum_Of_Weight` (Σw over the JEWEL events).
`Integrated_Lumi` = N/σ [pb⁻¹] follows the PHPythia8 convention. With `WEIGHTED T` a histogram filled
with w is normalised by **σ / Sum_Of_Weight**, where σ = N/L, as in `hwfilter`.

### Isospin (ISOCHANNEL XX)

Standalone JEWEL draws all four channel counts up front from NEVENT. A trigger-driven run does not
know how many events it will generate. PHJewel draws the channel per block of `set_isospin_block(n)`
JEWEL events (default 1000) with probabilities Z², Z(A−Z), (A−Z)Z, (A−Z)² over A². The mix is exact in
expectation, and PYTHIA is re-initialised only when the channel changes. σ is PARI(1) at the end of
each PYTHIA initialisation, weighted by the events made under it. For a fixed ISOCHANNEL the call
sequence is the standalone one, which is why seed 1000k reproduces NJOB k event for event.

### Truth trigger

The trigger is `HepMCTruthTriggerCore` from ProjectHERWIG: anti-k_T over status-1 particles without
neutrinos. With recoils it sees **unsubtracted** jets, because the status-3 thermal momenta are
excluded but the recoils are in. Those jets are harder than subtracted ones, so a trigger at X keeps
the events whose subtracted jets are above X.

## Files

| file | role |
|---|---|
| `PHJewel.{h,cc}` | the SubsysReco |
| `JewelCard.h` | strict card reader, canonical cards, hash |
| `JewelCore.h` | dlopen loader, event loop, isospin channels, bookkeeping (no Fun4All, no HepMC) |
| `JewelRun.h` | validate → canonical card → init → read-back check (shared by PHJewel and jwgen_test) |
| `jewel/jewel_capi.cc`, `jewel/jewel-api.f` | the core's C API and its Fortran side |
| `jewel/*.f` (generated), `pythia6425mod-lhapdf6.f`, `meix.f`, `subtraction.C` | JEWEL itself |
| `tools/make_api_sources.py` | generates the library sources from `../jewel-2.6.0` |
| `tools/compare_hepmc.py` | event-by-event comparison of two HepMC2 files |
| `tools/dst_to_hepmc.C` | DST → HepMC2 ASCII, plus PHGenIntegral and the first event's extras |
| `Fun4All_PHJewel.C` | macro |
| `jwgen_test.cc`, `Makefile.local` | laptop build and standalone test |

## For testing (Blair)

1. Same-seed regression on your machine. Run `bin/jewel-std-simple` on a card that has `NJOB k` and a
   fixed `ISOCHANNEL`. Run PHJewel with `set_seed(1000*k)` and `set_fortran_hepmc_file`. Then
   `cmp` the two HepMC files, and run `compare_hepmc.py` on `tools/dst_to_hepmc.C` output.
2. Geant: check that status-3 thermal partons and status-4 entries are **not** transported (they
   must not show up as G4 primaries) and that every status-1 particle is.
3. A recoil sample through the jet trigger, with subtraction at analysis level from the status-3
   particles in the DST truth record.
4. Normalisation: compare σ·Σ_sel w/Σw from a PHJewel DST with a standalone JEWEL histogram of the
   same card.
5. Try to break the cards. Every mistake JEWEL would swallow should be refused with a reason.
