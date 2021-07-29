#! /bin/env bash

readonly ARGS=("$@")
readonly ARGc="$#"
readonly ARGC=$(eval echo {1..$#..2})

readonly usage="

Usage:
Make manifest file [-h] [-d -o]
Arguments:
-h/--help shows this help text
-d|--reads_dir read directory
-o|--output output file name
"

error_message(){

	if [[ "$#" -eq 0  ]]
		then
			echo "Need to give an error message and usage"
		else
			printf "$1\n $2"
	fi

	exit 1
}

check_args(){
  echo "Running Argument Check"

	echo "${ARGS[0]}"


	#This is checking if any arguments were supplied
if [[ -z ${ARGS[0]} ]];
  then
  error_message "'Error: No arguments given.'  ${usage}"
fi
	#this is indexing through the -rst arguments and inputing them into the appropriate global variable.
	#arg_path is arg + 1 because it is the next argument.

	for arg in ${ARGC[@]};
		do
		arg_path=$(( ${arg}-1 ))
    # echo "${ARGS[$arg_path]}"
    # echo "${ARGS[$arg]}"
		case "${ARGS[$arg_path]}" in
		-h|--help) echo "${usuage}";;
		-d|--reads_dir) [[ -d "${ARGS[$arg]}" ]] && export dirpath=$(readlink -f ${ARGS[$arg]}) || error_message "'Error: Not a directory given' ${usage}";;
    -o|--output) [[ ! -z "${ARGS[$arg]}" ]] && export outname=${ARGS[$arg]} || "${outname}" || error_message "'Error: Not a file given' ${usage}";;
		*) printf "ERROR: incorrect argument uesd\n ${usage}";exit ;;
		esac

	done
}




create_manifest(){
  echo "man readpath: $dirpath"
  echo "man outpath: $outname"
  echo 'sample-id	forward-absolute-filepath	reverse-absolute-filepath' > $outname

  export readpath=$dirpath

  for forward in $(find -L $readpath -name "*R1*" -type f); do
      echo $forward
     reverse=${forward/R1/R2}

     samp_id_tmp=$(basename ${forward#*-})
     samp_id_tmp2=${samp_id_tmp%%_*}
     samp_id=${samp_id_tmp2%%-*}

     echo "$samp_id	$forward	$reverse" >> $outname

  done

}

main(){

  check_args
  echo "readpath: $dirpath"
  echo "outpath: $outname"
  create_manifest

}
main
