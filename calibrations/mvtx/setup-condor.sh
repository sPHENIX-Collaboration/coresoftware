#!/bin/bash

# first argument is install directory
if [ $# -ne 1 ]; then
  echo "Usage: $0 <install directory>"
  echo "defaulting to ${PWD}/install"
  mkdir -p ${PWD}/install
  INSTALLDIR=${PWD}/install
else
  INSTALLDIR=$1
fi

# if first argument is -h or --help, print help
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  echo "Usage: $0 <install directory>"
  exit 0
fi

echo "Installing MvtxCalibration in $1"
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh $INSTALLDIR

CURRENTDIR=$PWD
SRCDIR=$PWD/src
echo "MvtxCalibration installed in $INSTALLDIR"


# create condor scripts
echo "Creating condor scripts"
# create condor directory
mkdir -p $CURRENTDIR/condor
mkdir -p $CURRENTDIR/condor/log
cd $CURRENTDIR/condor
# see if condor job file exists
if [ -f "fhr_task.job" ]; then
  echo "fhr_task.job already exists"
  rm -f fhr_task.sh
fi 
echo "Creating fhr_task.job"
echo "Universe           = vanilla" > fhr_task.job
echo "initialDir         = ${CURRENTDIR}/condor" >> fhr_task.job
echo "Executable         = \$(initialDir)/run_fhr_task.sh" >> fhr_task.job
echo "PeriodicHold       = (NumJobStarts>=1 && JobStatus == 1)" >> fhr_task.job
echo "request_memory     = 6GB" >> fhr_task.job
echo "Priority           = 20" >> fhr_task.job
echo "Output             = \$(initialDir)/log/condor-mvtx-fhr-run$(run_number)-\$INT(Process,%05d).out" >> fhr_task.job
echo "Error              = \$(initialDir)/log/condor-mvtx-fhr-run$(run_number)-\$INT(Process,%05d).err" >> fhr_task.job
echo "Log                = /tmp/condor-mvtx-fhr-run$(run_number)-\$INT(Process,%05d)-\$ENV(USER)" >> fhr_task.job
echo "Arguments          = \$(nevents) \$(run_number) \$(trigger_rate) \$(input_file) \$(trigger_guard_file)" >> fhr_task.job
echo "Queue nevents run_number trigger_rate input_file trigger_guard_file from args.list" >> fhr_task.job

if [ -f "args.list" ]; then
    echo "args.list already exists"
    rm -f args.list
fi
echo "Creating args.list"
mkdir -p $CURRENTDIR/rootfiles
echo "1000000 42639 202 ${CURRENTDIR}/rootfiles/fhr_calib_42639_202kHz_1000000_events.root ${CURRENTDIR}/rootfiles/trigger_guard_42639_202kHz_1000000_events.root" > args.list
echo "1000000 42640 101 ${CURRENTDIR}/rootfiles/fhr_calib_42640_101kHz_1000000_events.root ${CURRENTDIR}/rootfiles/trigger_guard_42640_101kHz_1000000_events.root" >> args.list
echo "1000000 42641 44 ${CURRENTDIR}/rootfiles/fhr_calib_42641_44kHz_1000000_events.root ${CURRENTDIR}/rootfiles/trigger_guard_42641_44kHz_1000000_events.root" >> args.list
echo "-1 42639 202 ${CURRENTDIR}/rootfiles/fhr_calib_42639_202kHz_all_events.root none" >> args.list
echo "-1 42640 101 ${CURRENTDIR}/rootfiles/fhr_calib_42640_101kHz_all_events.root none" >> args.list
echo "-1 42641 44 ${CURRENTDIR}/rootfiles/fhr_calib_42641_44kHz_all_events.root none" >> args.list

# see if run_fhr_task.sh exists
if [ -f "run_fhr_task.sh" ]; then
    echo "run_fhr_task.sh already exists"
    rm -f run_fhr_task.sh
fi 

echo "#!/bin/bash" > run_fhr_task.sh
echo "" >> run_fhr_task.sh
echo "#arguments" >> run_fhr_task.sh
echo "if [ \$# -lt 1 ]; then" >> run_fhr_task.sh
echo "  echo \"Usage: \$0 <nevents> <runnumber> <trigger rate> <output file> <debug output file (optional)>\"" >> run_fhr_task.sh
echo "  exit 1" >> run_fhr_task.sh
echo "fi" >> run_fhr_task.sh
echo "" >> run_fhr_task.sh
echo "nevents=\$1" >> run_fhr_task.sh
echo "runnumber=\$2" >> run_fhr_task.sh
echo "trig_rate=\$3" >> run_fhr_task.sh
echo "outfile=\"\\\"\${4}\\\"\"" >> run_fhr_task.sh
echo "debugfile=\"\\\"\\\"\"" >> run_fhr_task.sh
echo "if [ \$5 == none ]; then" >> run_fhr_task.sh
echo "  debugfile=\"\\\"\\\"\"" >> run_fhr_task.sh
echo "else" >> run_fhr_task.sh
echo "  debugfile=\"\\\"\${5}\\\"\"" >> run_fhr_task.sh
echo "fi" >> run_fhr_task.sh
echo "" >> run_fhr_task.sh
echo "echo \"nevents: \${nevents}\"" >> run_fhr_task.sh
echo "echo \"runnumber: \${runnumber}\"" >> run_fhr_task.sh
echo "echo \"trig_rate: \${trig_rate}\"" >> run_fhr_task.sh
echo "echo \"outfile: \${outfile}\"" >> run_fhr_task.sh
echo "echo \"debugfile: \${debugfile}\"" >> run_fhr_task.sh
echo "" >> run_fhr_task.sh
echo "" >> run_fhr_task.sh
echo "source /opt/sphenix/core/bin/sphenix_setup.sh -n new" >> run_fhr_task.sh
echo "export INSTALLDIR=${INSTALLDIR}" >> run_fhr_task.sh
echo "source /opt/sphenix/core/bin/setup_local.sh \$INSTALLDIR" >> run_fhr_task.sh
echo "" >> run_fhr_task.sh
echo "cd ${CURRENTDIR}/macros" >> run_fhr_task.sh
echo "root -l -q -b Fun4All_MvtxFHR.C\(\${nevents},\${runnumber},\${trig_rate},\${outfile},\${debugfile}\)" >> run_fhr_task.sh
echo "" >> run_fhr_task.sh
echo "echo \$?" >> run_fhr_task.sh
echo "echo \"Script done\"" >> run_fhr_task.sh

# make run_fhr_task.sh executable
chmod +x run_fhr_task.sh

echo "Condor scripts created"
echo "Done"
cd $CURRENTDIR
