#!/bin/bash
##
## Functions for doing WWW Cross section measurement
## Bobak Hashemi

HIST_OUTPUT_LOCATION=`cat ConfigHelper.C | grep "^TString HIST_OUTPUT_LOCATION=" |sed "s/.*HIST_OUTPUT_LOCATION=\"\(.*\)\";*/\1/g"`
PLOT_OUTPUT_LOCATION=`cat ConfigHelper.C | grep "^TString PLOT_OUTPUT_LOCATION=" | sed "s/.*PLOT_OUTPUT_LOCATION=\"\(.*\)\";*/\1/g"`

function makePlots {
	if [[ $# < 1 ]]
	then
		echo "makePlots /path/to/config <no>(optional for no debug plots)"
		return
	fi
	mkdirs $1 plots
	conf_tmp_path=${1//.conf/.conf_tmp}
	./preprocessConf.py $1
	if [[ $2 == "no" ]]
	then
	    nice -n 19 root -l -b -q "drawPlots.C(\"$conf_tmp_path\", 0)" | tee ${PLOT_OUTPUT_LOCATION}${SR_IDENTITY}`basename $conf_filename .conf`/output.txt
	else
		nice -n 19 root -l -b -q "drawPlots.C(\"$conf_tmp_path\")" | tee ${PLOT_OUTPUT_LOCATION}${SR_IDENTITY}`basename $conf_filename .conf`/output.txt
	fi
	echo https://github.com/bth5032/WWWCrossSection/commit/`git rev-parse HEAD`/ > ${PLOT_OUTPUT_LOCATION}${SR_IDENTITY}`basename $conf_filename .conf`/commiturl.txt
}

function makeHistos {	
	mkdirs $2 hists
	conf_tmp_path=${2//.conf/.conf_tmp}
	./preprocessConf.py $2
	nice -n 19 root -l -b -q "doAll.C+(\"$1\", \"$conf_tmp_path\")" | tee ${HIST_OUTPUT_LOCATION}${SR_IDENTITY}/$1.output
	echo https://github.com/bth5032/WWWCrossSection/commit/`git rev-parse HEAD`/ > ${HIST_OUTPUT_LOCATION}${SR_IDENTITY}/commiturl.txt
}	

function setOutputLocations {
	if [[ $# > 0 ]]
	then
		if [[ ${1:(-4)} == "conf" ]]
		then
			SR_IDENTITY=${1#*configs/}
			SR_IDENTITY=`dirname $SR_IDENTITY`/
		else
			echo "You must give the location to a .conf file."
			return 0
		fi
	fi
}

function mkdirs {
	conf_filename=$1
	setOutputLocations $conf_filename 

	if [[ $2 == "hists" ]]
	then
		#Make Hist output location if it's not there
		if [[ ! -d ${HIST_OUTPUT_LOCATION}${SR_IDENTITY} ]]
		then	
			mkdir -p ${HIST_OUTPUT_LOCATION}${SR_IDENTITY}
		fi

		#Legacy Style For Hist Configs
		new_dir=`grep DEFAULT::histo_output_dir < $conf_filename`
		if [[ ! -z $new_dir ]]
		then
			mkdir -p ${new_dir#*=}
		fi
	fi

	if [[ $2 == "plots" ]]
	then
		#Make plot output location if it's not there.
		if [[ ! -d ${PLOT_OUTPUT_LOCATION}${SR_IDENTITY}`basename $conf_filename .conf`"/Debug/" ]]
		then
			mkdir -p ${PLOT_OUTPUT_LOCATION}${SR_IDENTITY}`basename $conf_filename .conf`"/Debug/" #Make up to path/to/config/Debug
			addIndexToDirTree ${PLOT_OUTPUT_LOCATION}${SR_IDENTITY}`basename $conf_filename .conf`"/Debug/" #Add index.php to each new folder
		fi

		#Legacy Style For Plots Configs
		new_dir=`grep DEFAULT::save_dir < $conf_filename`
		for l in `echo $new_dir`
		do
			if [[ ! -d ${l#*=} ]]
			then
				mkdir -p ${l#*=}"/Debug"
				addIndexToDirTree ${l#*=}"/Debug"
			fi
		done
	fi
}

function _makeAllForDir {
	fname_hist=${1//\//_} #turn path into _ seperated
	fname_hist=${fname_hist//__/_} #remove double underline string
	fname_hist=${fname_hist%_}.hist_out #remove trailing _, add extension

	fname_plots=${1//\//_} #turn path into _ seperated string
	fname_plots=${fname_plots//__/_} #remove double underline
	fname_plots=${fname_plots%_}.plots_out #remove trailing _, add extension

	if [[ $2 == "hists" ]]
	then
		echo $1 > outputs/$fname_hist
		makeHistosForDir $1 $3 >> outputs/$fname_hist 2>&1
	elif [[ $2 == "plots" ]]
	then
		echo $1 > outputs/$fname_plots
		makePlotsForDir $1 >> outputs/$fname_plots 2>&1
	elif [[ $2 == "all" ]]
	then
		echo $1 > outputs/$fname_hist
		makeHistosForDir $1 $3 >> outputs/$fname_hist 2>&1

		echo $1 > outputs/$fname_plots
		makePlotsForDir $1 >> outputs/$fname_plots 2>&1
	else
		echo "Please choose a step to run over, hists, plots, or all"
	fi
}

function makeAllForDir {
	if [[ $# < 2 ]]
	then
		echo "makeAllForDir <path_to_configs> <all/hists/plots> <sample_name | only if given hists or all>"
	elif [[ `find $1/* -maxdepth 0 -type d | wc -l` > 0 ]]
	then
		for i in `find $1/* -maxdepth 0 -type d`
		do
			echo "recursing into $i"
			makeAllForDir $i $2 $3
		done
	else
		echo -n `basename $1`" -- "
		_makeAllForDir $1 $2 $3 &
	fi
}

function makeHistosForDir {
	#Takes at most 2 args, the first is the name of the directory to run the run_modes over
	#the second is the name of the sample config. If no sample config is given, it runs over "all"

	if [[ $# < 2 ]]
	then
		makeHistosForDir_whichHist=all
	else
		makeHistosForDir_whichHist=$2
	fi

	if [[ -a $1/run_modes.conf ]]
	then
		makeHistos $makeHistosForDir_whichHist $1/run_modes.conf
	else
		echo "Can not find $1/run_modes.conf"
	fi
}

function makePlotsForDir {
	for i in `ls ${1}/*plots*.conf`
	do
		makePlots $i
	done
}

function addIndexToDirTree {
	#Adds the file at ~/public_html/WWWCrossSection/index.php into everything inside of the ~/public_html/ClosureTests/ directory for the directory given as $1.
	topdir=$1

	while [[ ${topdir%WWWCrossSection*} == "/home/users/bhashemi/public_html/" ]]
	do
		
		if [[ ! -a ${topdir}/index.php ]]
		then
			cp ~/public_html/index.php ${topdir}/index.php
		fi

		topdir=`dirname $topdir`
	
	done
}

function numjobs {
	psout=`ps aux | grep bhashemi | grep "configs/"` 

	if [[ $1 == "v" ]]
	then 
		echo "$psout" | grep "root -l"
	fi

	echo "scale=3; "`echo "$psout" | wc  -l`"*(1/2) - 1/2" | bc;
}

function pullHists {

	if [[ -z $HIST_OUTPUT_LOCATION || -z $SR_IDENTITY ]]
	then
		echo "please run setOutputLocations <path_to_run_modes.conf> first"
		return 1
	fi

	scp uaf:${HIST_OUTPUT_LOCATION}${SR_IDENTITY}*.root histos/ 
}

function killjobs {
	kill -9 `ps aux | grep "^bhashemi" | grep "root" | head -n-1 | cut -d' ' -f2 | xargs`
}

function addHists {
	if [[ $# < 2 ]]
	then
		echo "Calls hadd for you and determines which files you mean and where to put output by file's config Name"
		echo "addHists <output_fname> <input_1> <input_2> ..."
		return 1
	fi

	if [[ -z $HIST_OUTPUT_LOCATION || -z $SR_IDENTITY ]]
	then
		echo "please run setOutputLocations <path_to_run_modes.conf> first"
		return 1
	fi

	HADD_LIST=${HIST_OUTPUT_LOCATION}${SR_IDENTITY}${2}".root"

	for arg in ${@:3}
	do
		HADD_LIST="$HADD_LIST ${HIST_OUTPUT_LOCATION}${SR_IDENTITY}${arg}.root" 
	done

	echo "Running: hadd ${HIST_OUTPUT_LOCATION}${SR_IDENTITY}$1.root $HADD_LIST"

	hadd ${HIST_OUTPUT_LOCATION}${SR_IDENTITY}$1.root $HADD_LIST
}

function backupConfig {
	#===========================================================
	# Backs up configuration with the given config directory, allows
	# you to leave a message in the backup dir which tells what has
	# changed in the new version.
	#===========================================================
	
	backupConfig_path=$1

	if [[ $# < 1 ]]
	then
		echo "backupConfig configs/path/to/confdir"
  	read -e -p "Config Dir Path: " backupConfig_path
  fi

  backupConfig_path=${backupConfig_path%/} #Remove Trailing '/' if it's there.

  SR_IDENTITY=${backupConfig_path#*configs/}

	backupConfig_movePlots="n"
	backupConfig_moveHists="n"
	backupConfig_timestamp=`date +%m%d%y_%Hh%M`
	backupConfig_message=

	backupConfig_histdir=${HIST_OUTPUT_LOCATION}${SR_IDENTITY} #Current hists dir
	backupConfig_histbak=${backupConfig_histdir}_bak${backupConfig_timestamp} #New hists dir

	backupConfig_plotdir=${PLOT_OUTPUT_LOCATION}${SR_IDENTITY} #Current plots dir
	backupConfig_plotbak=${backupConfig_plotdir}_bak${backupConfig_timestamp} #New plots dir

  if [[ ! -d $backupConfig_path ]]
	then
		echo "Could not find directory: "$backupConfig_path", please try again"
  	return
  fi

 	if [[ -d ${backupConfig_histdir} ]] 
 	then
 		read -p "move hists at ${backupConfig_histdir} to ${backupConfig_histbak}? (y/n)" backupConfig_moveHists
	else
		echo "No previous hists found at ${backupConfig_histdir}"
	fi

	if [[ -d ${backupConfig_plotdir} ]] 
 	then
 		read -p "move plots at ${backupConfig_plotdir} to ${backupConfig_plotbak}? (y/n)" backupConfig_movePlots
	else
		echo "No previous plots found at ${backupConfig_plotdir}"
	fi

	read -e -p "What changes in the next version: " backupConfig_message

	if [[ "$backupConfig_moveHists" == "y" ]]
	then
		cp -r $backupConfig_histdir $backupConfig_histbak
		echo $backupConfig_message > ${backupConfig_histbak}/backup_log.txt
		echo "Moved Hists"
	fi

	if [[ "$backupConfig_movePlots" == "y" ]]
	then
		cp -r $backupConfig_plotdir $backupConfig_plotbak
		echo $backupConfig_message > ${backupConfig_plotbak}/backup_log.txt
		echo "Moved Plots"
	fi
}
