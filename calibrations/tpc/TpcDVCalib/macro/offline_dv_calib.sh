path=$PWD

runnumber=$1
cd Reconstructed/$runnumber/

FILES=""
MERGE="final_${runnumber}_ana.root"
for f in TrackCalo_*_ana.root
do
  echo "Processing $f"
  FILES+="$f "
done
hadd -f -k $MERGE $FILES

cd $path

root -b -q -l saveroot.C\($runnumber\)

mkdir figure
root -b -q -l fit.C\($runnumber\)
