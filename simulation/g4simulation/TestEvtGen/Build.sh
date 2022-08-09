echo "-------------- CHANGING SOFTWARE NOW BRO --------------------"

#source /opt/sphenix/core/bin/sphenix_setup.csh -n ana.270
#source /cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.csh -n ana.221
#source /opt/sphenix/core/bin/sphenix_setup.csh -n ana.200

#source /opt/sphenix/core/bin/sphenix_setup.csh -n 

#source /cvmfs/sphenix.sdcc.bnl.gov/gcc-8.3/opt/sphenix/core/bin/sphenix_setup.csh -n ana.269


#source /opt/sphenix/core/bin/sphenix_setup.csh -n ana.269
#source /opt/sphenix/core/bin/sphenix_setup.csh -n 

rm -r install

mkdir install
setenv MYINSTALL $PWD/install/
setenv LD_LIBRARY_PATH $MYINSTALL/lib:$LD_LIBRARY_PATH
set path = ( $MYINSTALL/bin $path )



echo "-------------- BUILD main ------------------------"


#cd coresoftware/simulation/g4simulation/g4detectors/
sh autogen.sh --prefix=$MYINSTALL
make
make install
