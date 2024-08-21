path=$PWD

runnumber=$1

dir_name="figure"
if [ ! -d "$dir_name" ]; then
  mkdir -p "$dir_name"
  echo "Directory $dir_name created."
else
  echo "Directory $dir_name already exists."
fi

dir_name="figure_failed"
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

dir_name="cdbttree_failed"
if [ ! -d "$dir_name" ]; then
  mkdir -p "$dir_name"
  echo "Directory $dir_name created."
else
  echo "Directory $dir_name already exists."
fi

root -b -q -l fit.C\($runnumber\)
