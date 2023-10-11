#!/bin/bash

# delete the files to output values to
start=1 			# set to 1 in order to start from the beginning
num_threads=32 		# check 8, 16, 32
if [[ $start == 1 ]]; then
	rm -f outfile.log
	rm -f outfile_raw.log
	rm -f outfile_joined.log
fi
echo $num_threads threads >> outfile.log
echo $num_threads threads >> outfile_raw.log

# run the tests
for str in $(ls soplib);
do
	if [[ $start == 1 ]]; then
		if [[ $str != "R.400.1000.1.sop" && $str != "R.600.1000.1.sop" && $str != "R.600.1000.15.sop" && $str != "R.700.1000.1.sop" && $str != "R.700.1000.15.sop" ]]; then
			echo $str >> outfile.log
			output=$(timeout -v 1h ./sop_solver soplib/$str 8 soplib_config.txt)
			echo "$str" >> outfile_raw.log
			echo "$output" >> outfile_raw.log
			exit_status=$?
			if [[ $exit_status -eq 124 ]]; then # test times out
				echo TIMED OUT >> outfile.log
			elif [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # test produces no usable timing data
					echo NO TIMING DATA >> outfile.log
			else # output timing data to file
				echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile.log # add " | cut -d, -f2" to find only the time
			fi
		fi
	elif [[ $str == "R.300.100.60.sop" ]]; then
		echo "Starting..."
		start=1
	fi 
done

# join pairs of lines for processing by gnuplot
cat outfile.log | paste -d " " - - >> outfile_joined.log
