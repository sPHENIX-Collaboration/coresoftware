#! /bin/bash
verbose_mode=false
events=1000000
nfiles=1000
dosubmit=false
triggertype="MB" 
configfile="MB.in"
configdir="$(pwd)/../config_files"
condor_testfile="condor_blank.job"
user=`id -u -n`
make_condor_jobs()
{
	for i in $(seq 0 ${nfiles}); do 
		condor_file="condor_file_dir/condor_"$triggertype"_"$i".job"
		condor_out_file="condor_file_dir/condor_"$triggertype"_"$i".out"
		condor_err_file="condor_file_dir/condor_"$triggertype"_"$i".err"
		condor_log_file="condor_file_dir/condor_"$triggertype"_"$i".log"
		if [ "$vebose_mode" = true ]; then
			echo "Producing condor job file " $condor_file
		fi
		IFS=$'\n' read -d '' -r -a blanklines < $condor_testfile
		echo "${blanklines[0]}" > $condor_file 
		echo "${blanklines[1]}"$(pwd)"/Herwig_run.sh" >> $condor_file
		echo "${blanklines[2]}"$configfile "  1000 " $i >> $condor_file
		echo "${blanklines[3]}"$condor_out_file >> $condor_file
		echo "${blanklines[4]}"$condor_err_file >> $condor_file
		echo "${blanklines[5]}"$condor_log_file >> $condor_file
		echo "${blanklines[6]}"$(pwd) >>$condor_file
		echo "${blanklines[7]}" >> $condor_file
		echo "${blanklines[8]}" >> $condor_file #set for me right now
		echo "${blanklines[9]}" "   "  $user >> $condor_file 
		echo "${blanklines[10]}" >> $condor_file
		echo "${blanklines[11]}" >> $condor_file
		echo "${blanklines[12]}" >> $condor_file
		echo "${blanklines[13]}" >> $condor_file
	done		
}
submit_condor_jobs(){
	#if submit just get all files in the expected job type	
	ls condor_file_dir
	for i in `ls "condor_file_dir/condor_"$triggertype"_"*".job"`; do 
		condor_submit $i 
	done
}
has_argument(){
	[[ ("$1" == *=* && -n ${1#*=}) || ( ! -z "$2" && "$2" != -*) ]]
}

extract_argument() {
	echo "${2:-${1#*=}}"
}

set_config()
{
	if [ "$triggertype" = "MB" ]; then
		configfile="${configdir}/MB.in"
	elif [ "$triggertype" = "Jet10" ]; then
		configfile="${configdir}/Jet10.in"
	elif [ "$triggertype" = "Jet30" ]; then 
		configfile="${configdir}/Jet30.in"
	else
		configfile="${configdir}/MB.in" #use as default value
	fi
}

handle_options(){
	while [ $# -gt 0 ]; do
	case $1 in 
		-h | --help)
			echo "Options for Herwig job creation script"
			echo "$0 [OPTIONS]"
			echo "This script runs Herwig to create HepMC files given an input configuration"
			echo " \n "
			echo " -h, --help	Display this help message"
			echo " -v, --verbose	Enable verbose job creation"
			echo " -N, --events  	Number of events to generate"
			echo " -s, --submit	Make and submit condor jobs"
			echo " -t, --trigger	Input type (MB, Jet10, Jet30)"
			echo " -i, --input 	Specify new input file"
			exit 0 
			;;
		-v | --verbose)
			verbose_mode=true
			shift
			;;
		-N | --events*)
			if has_argument $@; then 
				events=$(extract_argument $@) 
				nfiles=$(( events / 1000 ))
				if [ "$verbose_mode" = true ]; then
					echo "Run " $events " events"
					echo " This will generate " $nfiles " output hepmc files"
				fi
				nfiles=$(( nfiles - 1 ))
			fi
			shift
			shift
			;;
		-s | --submit) 
			dosubmit=true
			shift			
			;;
		-t | --trigger*)
			if has_argument $@; then 
				triggertype=$(extract_argument $@)
				set_config
				if [ "$verbose_mode" = true ]; then
					echo "Trigger type: " $triggertype 
					echo "Config file: " $configfile
					echo "Config dir: " $configdir
				fi
			fi
			shift
			shift
			;;
		-i | --input*)
			if has_argument $@; then
				configfile=$(extract_argument $@)
				if [ "$verbose_mode" = true]; then
					echo "Trigger type: " $triggertype 
					echo "Config file: " $configfile
				fi
			fi
			shift
			;;
		*) 
			echo "Invalid option: $1 "
			exit 1
			;;
		esac
	done
}

handle_options "$@"
make_condor_jobs 
if [ "$dosubmit" = true ]; then
	submit_condor_jobs
fi

