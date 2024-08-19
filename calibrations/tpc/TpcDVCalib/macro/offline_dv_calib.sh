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

dir_name="figure"
if [ ! -d "$dir_name" ]; then
  mkdir -p "$dir_name"
  echo "Directory $dir_name created."
else
  echo "Directory $dir_name already exists."
fi

dir_name="root"
if [ ! -d "$dir_name" ]; then
  mkdir -p "$dir_name"
  echo "Directory $dir_name created."
else
  echo "Directory $dir_name already exists."
fi

dir_name="cdbttree"
if [ ! -d "$dir_name" ]; then
  mkdir -p "$dir_name"
  echo "Directory $dir_name created."
else
  echo "Directory $dir_name already exists."
fi

root -b -q -l fit.C\($runnumber\)
