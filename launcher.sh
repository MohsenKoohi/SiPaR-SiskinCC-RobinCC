#!/bin/bash
#
# This script runs an algorithm for a number of graph datasets


function __expand_hostlist() {
	local hostlist=$1
	if [[ $hostlist != *\[* ]]; then
		echo -n "$hostlist,"
		return
	fi
	local prefix=${hostlist%%[*}
	local ranges=${hostlist#*[}
	ranges=${ranges%]*}

	for range in ${ranges//,/ }; do
		if [[ $range =~ - ]]; then
			local start=${range%-*}
			local end=${range#*-}
			for ((i=start; i<=end; i++)); do
				echo -n "${prefix}${i},"
			done
		else
			echo -n "${prefix}${range},"
		fi
	done
}

function hosts_list() {
	local output=""
	for input in "$@"; do
		output+=$(__expand_hostlist $input | sed 's/,$//')","
	done
	echo ${output%,}
}

function __count_hosts() {
	local hostlist=$1
	local prefix=${hostlist%%[*}
	local ranges=${hostlist#*[}
	ranges=${ranges%]*}
	local count=0

	for range in ${ranges//,/ }; do
		if [[ $range =~ - ]]; then
			local start=${range%-*}
			local end=${range#*-}
			count=$((count + end - start + 1))
		else
			count=$((count + 1))
		fi
	done
	echo $count
}

function hosts_count() {
	local total=0
	for input in "$@"; do
		total=$((total + $(__count_hosts $input)))
	done
	echo $total
}

function ul2s()
{
	in=$1
	ul2s_out=""

	if [ `echo "$in > 10^21" | bc` = 1 ]; then 
		ul2s_out=`echo "scale=3; $in/10^21" | bc`"Z"
	else if [ `echo "$in > 10^18" | bc` = 1 ]; then 
		ul2s_out=`echo "scale=3; $in/10^18" | bc`"E"
	else if [ `echo "$in > 10^15" | bc` == 1 ]; then 
		ul2s_out=`echo "scale=3; $in/10^15" | bc`"P"
	else if [ `echo "$in > 10^12" | bc` == 1 ]; then 
		ul2s_out=`echo "scale=3; $in/10^12" | bc`"T"
	else if [ `echo "$in > 10^9" | bc` == 1 ]; then 
		ul2s_out=`echo "scale=3; $in/10^9" | bc`"G"
	else if [ `echo "$in > 10^6" | bc` == 1 ]; then 
		ul2s_out=`echo "scale=3; $in/10^6" | bc`"M"
	else if [ `echo "$in > 10^3" | bc` == 1 ]; then 
		ul2s_out=`echo "scale=3; $in/10^3" | bc`"k"
	else 
		ul2s_out=$in
	fi fi fi fi fi fi fi

	return;
}

# Print help
	echo -e "\n\033[1;35mDistributed LaganLighter Launcher Script\033[0;37m\n"
	echo "Arguments to be passed to this script:"
	echo "  alg=algorithm filename (with or without .c)"
	echo "  df=path/to/datasets/folder (default ./data)"
	echo "  make-flags=\"flags (e.g. \`threads\`, \`no_ht\`, \`energy\`, \`papi\`,\`debug\`, and/or \`wait_active\`) to be passed to the make, separated by space, and double quoted\" (default empty)"
	echo "  program-args=\"arguemnts to be passed to the \`program\` by \`make\`, separated by space, and double quoted\" (default empty)"
	echo "  sf=m (start from dataset m, numbers starting from 0, default: 0)"
	echo "  sa=n (stop after processing dataset n, numbers starting from 0, default: -1)"
	echo "  -ld (just list datasets, default: 0)"
	echo "  -iw (include weighted graphs, default: false)"
	echo "  -shm-store (store graphs in shm, default: 0)"
	echo "  -shm-delete (delete shm graphs at the end, default: 0)"
	echo "  report-title=\"a double quoted string that will be printed on top of the report file\" (default empty)"
	echo "  mpiexec-args=\"arguments to mpiexec, double quoted\" (default \" --report-bindings -n \$SLURM_NNODES -npernode 1 -oversubscribe\") "

	echo

# Reading hosts
	hosts_packed=`scontrol show job $SLURM_JOB_ID | grep -e "\bNodeList" | cut -d '=' -f 2- | xargs`
	nodes_count=`hosts_count $hosts_packed`
	nodes_list=`hosts_list $hosts_packed`
	echo "Nodes count: "$nodes_count
	echo "Nodes list: "$nodes_list

# Reading arguments
	echo -e "\033[0;34mArguments\033[0;37m"
	ALG=`echo $@ | grep -oP "(alg=?[\S]*)"`
	if [ `echo $ALG | grep -oc "="` -eq 1 ]; then
		ALG=`echo $ALG | cut -f2 -d=`;
	fi
	if [ `echo $ALG | grep -oc "\.c"` -eq 1 ]; then
		ALG=`echo $ALG | cut -f1 -d.`;
	fi
	if [ ! -e "$ALG.c" ]; then 
		echo -e "  \033[0;31mError:\033[0;37m File \"$ALG.c\" does not exist.";
		exit 1;
	fi
	echo "  Algorithm: $ALG"

	DF="./data"
	MAKE_FLAGS=""
	PROGRAM_ARGS=""
	IW=0
	SF=0
	SA="-1"
	LD=0
	SHM_STORE=0
	SHM_DELETE=0
	REPORT_TITLE=""
	MPIEXEC_ARGS="--report-bindings -n $nodes_count -npernode 1 -oversubscribe" # --host $nodes_list

	for i in `seq 1 $#`; do
		if [[ "${!i}" == *"df"* ]]; then 
			DF=`echo "${!i}" | cut -f2 -d=`
		fi

		if [[ "${!i}" == *"make-flags"* ]]; then 
			MAKE_FLAGS=`echo "${!i}" | cut -f2- -d=`
		fi

		if [[ "${!i}" == *"program-args"* ]]; then 
			PROGRAM_ARGS=`echo "${!i}" | cut -f2- -d=`
		fi

		if [[ "${!i}" == *"sf="* ]]; then 
			SF=`echo "${!i}" | cut -f2 -d=`
		fi

		if [[ "${!i}" == *"sa="* ]]; then 
			SA=`echo "${!i}" | cut -f2 -d=`
		fi

		if [ "${!i}" == "-ld" ]; then
			LD=1
		fi

		if [ "${!i}" == "-iw" ]; then
			IW=1
		fi

		if [ "${!i}" == "-shm-store" ]; then
			SHM_STORE=1
		fi

		if [ "${!i}" == "-shm-delete" ]; then
			SHM_DELETE=1
		fi

		if [[ "${!i}" == *"report-title"* ]]; then 
			REPORT_TITLE=`echo "${!i}" | cut -f2- -d=`
		fi

		if [[ "${!i}" == *"mpiexec-args"* ]]; then 
			MPIEXEC_ARGS=`echo "${!i}" | cut -f2- -d=`
		fi

	done 

	if [ ! -e $DF ]; 
	then
		echo -e "  \033[0;31mError:\033[0;37m Datasets folder \"$DF\" does not exist.";
		exit 2;
	fi;

	
	
	echo "  Make flags (make-flags): \"$MAKE_FLAGS\""
	echo "  Program args (program-args): \"$PROGRAM_ARGS\""
	echo "  Datasets folder (df): $DF" 
	echo "  Start from dataset (sf): $SF"
	echo "  Stop after dataset (sa): $SA"
	echo "  List datasets (-ld): $LD"
	echo "  Include weighted graphs (-iw): $IW"
	echo "  Store graph in shm (-shm-store): $SHM_STORE"
	echo "  Delete shm graphs at end (-shm-delete): $SHM_DELETE"
	echo "  Report title (report-title): $REPORT_TITLE"
	echo "  mpiexec args (mpiexec-args): $MPIEXEC_ARGS"
	echo 

# Identifying datasets
	echo -e "\033[0;34mDatasets\033[0;37m"
	DF=`realpath $DF`
	echo "  Dataset folder realpath: $DF"
	dataset_patterns="*.graph *.txt"
	if [ $IW == 1 ]; then
		dataset_patterns="$dataset_patterns *.labels"
	fi

	initial_datasets=(`(cd $DF; ls $dataset_patterns 2>/dev/null )`)
	datasets=""
	dataset_edges=""
	dataset_vertices=""
	for (( i=0; i <  ${#initial_datasets[@]}; i++)); do
		ds=${initial_datasets[$i]}
		
		suffix=`echo $ds | rev | cut -f1 -d. | rev`
		pre_suffix=`echo $ds | rev | cut -f1 -d_ | rev`

		edges_count=0
		vertices_count=0

		if [ $pre_suffix == "props.txt" ]; then
			[ ! -e $DF/`echo $ds | sed 's/props.txt/offsets.bin/'` ] && continue
			[ ! -e $DF/`echo $ds | sed 's/props.txt/edges.bin/'` ] && continue

			vertices_count=`cat "$DF/$ds" | grep vertices-count | cut -f2 -d:`
			edges_count=`cat "$DF/$ds" | grep edges-count | cut -f2 -d:`
			ds=`echo $ds | sed 's/_props.txt/.bin/g'`
		elif [ $suffix == "txt" ]; then
			edges_count=`head "$DF/$ds" -n2 | tail -n1`
			vertices_count=`head "$DF/$ds" -n1`
		fi

		if [ $suffix == "graph" ]; then
			prop_file=`echo $ds | sed 's/\(.*\)\.graph/\1\.properties/'`
			[ ! -e $DF/$prop_file ] && continue
			edges_count=`cat "$DF/$prop_file" | grep -P "\barcs[\s]*=" | cut -f2 -d= | xargs`
			vertices_count=`cat "$DF/$prop_file" | grep -P "\bnodes[\s]*=" | cut -f2 -d= | xargs`
		fi

		if [ $suffix == "labels" ]; then
			laebls_prop_file=`echo $ds | sed 's/\(.*\)\.labels/\1\.properties/'`
			prop_file=`cat "$DF/$laebls_prop_file" | grep -P "underlyinggraph[\s]*=" | cut -f2 -d=|xargs`".properties"
			[ ! -e $DF/$prop_file ] && continue
			edges_count=`cat "$DF/$prop_file" | grep -P "\\barcs=" | cut -f2 -d=`
			edges_count=`cat "$DF/$prop_file" | grep -P "\\bnodes=" | cut -f2 -d=`
		fi

		if [ `echo $edges_count |  tr -cd ' \t' | wc -c ` != 0 ]; then
			edges_count=0;
		fi

		[[ ! -n "$edges_count" || $edges_count == "0" ]] &&
			echo -e "  \033[0;31mError\033[0;37m: dataset \"$DF/$ds\" could not be recognized." &&
			continue

		datasets="$datasets $ds"
		dataset_edges="$dataset_edges $edges_count"
		dataset_vertices="$dataset_vertices $vertices_count"

		# echo $i $ds $suffix $edges_count $vertices_count;
	done	

	datasets=($datasets);
	dataset_edges=($dataset_edges);
	dataset_vertices=($dataset_vertices);

	sorted_datasets=`for (( i=0; i < ${#datasets[@]}; i++)); do
		echo ${datasets[$i]} ${dataset_vertices[$i]} ${dataset_edges[$i]}
	done | sort -n -k3`

	datasets=""
	c=0
	while IFS= read -r line; do
		if [ $c -lt $SF ]; then
			c=`echo "$c+1"|bc`
			continue
		fi

		ds=`echo $line|cut -f1 -d' '`
		vc=`echo $line|cut -f2 -d' '`
		ec=`echo $line|cut -f3 -d' '`
		ul2s $vc
		vc_out=$ul2s_out
		ul2s $ec
		ec_out=$ul2s_out
		printf "  %2d- %-30s, |V|: %10s, |E|: %10s\n" $c $ds $vc_out $ul2s_out
		unset vc_out ec_out ul2s_out 
		datasets="$datasets $ds"

		if [ $SA -eq $c ]; then
			break
		fi
		c=`echo "$c+1"|bc`
		
	done <<< "$sorted_datasets"

	if [ $LD -eq "1" ]; then
		echo -e "\n"
		exit
	fi

# Dropping cache
	echo -e "\n\033[0;34mDropping OS cache\033[0;37m"
	if [[ $SA != 0 ]]; then
		env mpiexec $MPIEXEC_ARGS make clean 1>/dev/null 2>&1
	fi;
	env make 1>/dev/null 2>&1

# Processing graphs
	echo -e "\n\033[0;34mProcessing Datasets\033[0;37m"

	# Preparing the build folder
	hostnames=`echo $hosts_packed | sed -r 's/\[|\]//g' | sed -r 's/ /,/g'`
	home_folder=`pwd`
	build_folder=$home_folder/build/`date +"%Y%m%d-%H%M%S"`-$hostnames
	echo "  Build folder: $build_folder"
	mkdir -p $build_folder
	
	for f in `ls`; do
		[ $f = "build" ] && continue
		[ $f = "data" ] && continue
		[ $f = "obj" ] && continue
		[ $f = "logs" ] && continue
		[ $f = "docs" ] && continue
		[ $f = "others" ] && continue
		# [ $f = "paragrapher" ] && continue

		cp -r $f $build_folder
	done

	if [ ! -d logs ]; then
		mkdir logs
		c=0
	else
		c=`(cd logs; find  -maxdepth 1 \( -name "$ALG-$hostnames*" \) | wc -l)`
	fi
	log_folder="logs/$ALG-$hostnames-$c"

	while [ -d $log_folder ]; do
		c=`echo "$c+1"|bc`
		log_folder="logs/$ALG-$hostnames-$c"
	done

	mkdir -p $log_folder
	echo "  Log Folder: "$log_folder
	
	report_path=$log_folder/report.txt
	echo "  Report file: "$report_path

	# Header of report file
		echo "$ALG, $REPORT_TITLE" > $report_path
		echo -e "Machines: $hostnames\nGCC: `gcc -dumpversion`" >> $report_path
		echo -e "mpiexec: `mpiexec --version | xargs`" >> $report_path
		echo -e "mpiexec-args: $MPIEXEC_ARGS" >> $report_path
		echo -e "make-flags: "$MAKE_FLAGS"\nprogram-args: "$PROGRAM_ARGS >> $report_path
		echo "Datasets: $SF - $SA" >> $report_path
		c=$SF;
		for ds in $datasets; do
			echo "  $c, $ds" >> $report_path
			c=`echo "$c+1"|bc`
		done
		echo -e "\nStarted on: "`date +"%Y/%m/%d-%H:%M:%S"`"\n" >> $report_path
		echo

	c=$SF;
	for ds in $datasets; do

		suffix=`echo $ds | rev | cut -f1 -d. | rev`
		input_graph=""
		input_type=""

		if [ $suffix == "bin" ]; then
			input_graph="$DF/"`echo $ds | sed 's/\(.*\)\.bin/\1/'`
			input_type="binary"
		elif [ $suffix == "txt" ]; then
			input_graph="$DF/$ds"
			input_type="text"
		elif [ $suffix == "graph" ]; then
			input_graph="$DF/"`echo $ds | sed 's/\(.*\)\.graph/\1/'`
			input_type="PARAGRAPHER_CSX_WG_800_AP"
		elif [ $suffix == "labels" ]; then
			input_graph="$DF/"`echo $ds | sed 's/\(.*\)\.labels/\1/'`
			input_type="PARAGRAPHER_CSX_WG_404_AP"
		fi

		is_sym=0
		if [[ $ds == "sym"* ]]; then
			is_sym=1
		fi
		if [ `echo $ds | grep -ocP "MS[\d]+"` == "1"  ]; then
			is_sym=1
		fi

		ds_log_file="$log_folder/$c-"`echo $ds | rev | cut -f2- -d. | rev`".txt"

		cmd="LL_INPUT_GRAPH_PATH=$input_graph LL_INPUT_GRAPH_TYPE=$input_type"
		cmd="$cmd LL_INPUT_GRAPH_BATCH_ORDER=$c LL_INPUT_GRAPH_IS_SYMMETRIC=$is_sym"
		cmd="$cmd LL_STORE_INPUT_GRAPH_IN_SHM=$SHM_STORE LL_OUTPUT_REPORT_PATH=$home_folder"/"$report_path"
		cmd="$cmd LL_GRAPH_RA_BIN_FOLDER=$DF/RA_bin_arrays"
		
		echo "  $c, $input_graph:"
		echo "    $cmd make $ALG $MAKE_FLAGS mpiexec-args=\"$MPIEXEC_ARGS\" args=\"$PROGRAM_ARGS\" 1>$ds_log_file 2>&1" 

		cd $build_folder
		env $cmd make $ALG $MAKE_FLAGS mpiexec-args="$MPIEXEC_ARGS" args="$PROGRAM_ARGS" 1>$home_folder"/"$ds_log_file 2>&1
		cd $home_folder

		echo 
		
		c=`echo "$c+1"|bc`
	done

	rm -rf $build_folder

	echo -e "\nFinished on: "`date +"%Y/%m/%d-%H:%M:%S"`"\n" >> $report_path

# Post processing actions 
	if [ $SHM_DELETE == 1 ]; then 
		make clean_shm_graphs
	fi

# Printing the report
	echo -e "\n\033[0;34mResults\033[0;37m"
	echo "  Report file: $report_path"
	echo 
	while IFS= read -r line; do
		echo "  $line" 
	done <<< `cat $report_path`

echo 