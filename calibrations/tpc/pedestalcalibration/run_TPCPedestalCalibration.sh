source /opt/sphenix/core/bin/sphenix_setup.sh -n new

inputFiles="{"
for fileList in $@
do
  inputFiles+="\"${fileList}\","
done
inputFiles=${inputFiles::-1}
inputFiles+="}"
echo running: run_TPCPedestalCalibration.sh $*
root.exe -q -b TPCPedestalCalibration.C\(${inputFiles}\)
echo Script done
