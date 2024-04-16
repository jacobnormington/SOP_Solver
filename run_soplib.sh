#!/bin/bash

#if startInstance is set, it will add to the existing log files, picking up from a previous pause
startInstance=""					# the first instance to process, skipping any lexicographically earlier; if empty, begins at the first instance
endInstance=""						# the last instance to process, skipping any lexicographically later; if empty, continues through the entire suite
dataset=1							# which instances to consider: 1 to exclude very hard instances, 0 for all instances, -1 for only very hard instances
num_threads=32 						# check 8, 16, 32



start=0 			# DO NOT USE, use startInstance instead
config="soplib_config.txt"
if [[ $startInstance == "" ]]; then
	start=1
	rm -f outfile_sop.log
	rm -f outfile_sop_raw.log
	rm -f outfile_sop_joined.log

	echo $num_threads threads >> outfile_sop.log
	echo $num_threads threads >> outfile_sop_raw.log
fi



# run the tests
for str in $(ls -v soplib);
do
	if [[ $str == $startInstance ]]; then
		echo "Starting..."
		start=1
	fi

	if [[ $start == 1 ]]; then
		if [[ ${str: -4:4} == ".sop" ]]; then #discard anything other than sop instance files
			if [[ 	$dataset == 0 ||
					($dataset -gt 0 && ($str != "R.400.1000.1.sop" && $str != "R.600.1000.1.sop" && $str != "R.600.1000.15.sop" && $str != "R.700.1000.1.sop"
						&& $str != "R.700.1000.15.sop")) ||
					($dataset -lt 0 && !($str != "R.400.1000.1.sop" && $str != "R.600.1000.1.sop" && $str != "R.600.1000.15.sop" && $str != "R.700.1000.1.sop" 
						&& $str != "R.700.1000.15.sop"))
				]]; then
				echo $str >> outfile_sop.log
				output=$(./sop_solver soplib/$str $num_threads $config)
				echo "$str" >> outfile_sop_raw.log
				
				echo "$output" >> outfile_sop_raw.log

				if [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # time and cost not printed at the end
					echo NO TIMING DATA >> outfile_sop.log
				else # output timing data to file
					echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile_sop.log # add " | cut -d, -f2" to find only the time
				fi

				if [[ $( echo "$output" | egrep -o "instance timed out" ) != "" ]]; then # test times out
					echo TIMED OUT >> outfile_sop.log
				fi

				# if [[ $( echo $output | egrep -o "Total Progress = 100%" ) == "" ]]; then # if test doesn't end in 100% progress
				# 	echo PROGRESS BAR NOT 100% >> outfile_sop.log
				# 	echo PROGRESS BAR NOT 100% >> outfile_sop_raw.log
				# fi
			fi
		fi
	fi

	if [[ $str == $endInstance ]]; then
		echo "...script ended early"
		exit 1
	fi
done

# join pairs of lines for processing by gnuplot
# cat outfile_sop.log | paste -d " " - - >> outfile_joined.log
