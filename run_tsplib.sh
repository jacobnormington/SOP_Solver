#!/bin/bash

# delete the files to output values to
start=1 			# set to 1 in order to start from the beginning
num_threads=32 		# check 8, 16, 32
if [[ $start == 1 ]]; then
	rm -f outfile_tsp.log
	rm -f outfile_tsp_raw.log
	rm -f outfile_joined_tsp.log
fi
echo $num_threads threads >> outfile_tsp.log
echo $num_threads threads >> outfile_tsp_raw.log

# run the tests
for str in $(ls tsplib);
do
	if [[ $start == 1 ]]; then
		if [[ $str != "ESC07_last.sop" && $str != "ESC11_last.sop" ]]; then #these are not "real" instances, ask Taspon about them
			if [[ $str != "ESC78.sop" && $str != "ft53.1.sop" && $str != "ft53.2.sop" && $str != "ft53.3.sop" && $str != "ft70.1.sop" && $str != "ft70.2.sop" && $str != "ft70.3.sop" 
				&& $str != "kro124p.1.sop" && $str != "kro124p.2.sop" && $str != "kro124p.3.sop" && $str != "kro124p.4.sop" && $str != "p43.1.sop" && $str != "p43.2.sop" 
				&& $str != "prob.100.sop" && $str != "rbg323a.sop" && $str != "rbg341a.sop" && $str != "rbg358a.sop" && $str != "rbg378a.sop" 
				&& $str != "ry48p.1.sop" && $str != "ry48p.2.sop" && $str != "ry48p.3.sop" ]]; then #skip the instances that always timeout on 32 threads
				echo $str >> outfile_tsp.log
				output=$(timeout -v 1h ./sop_solver tsplib/$str $num_threads tsplib_config.txt)
				echo "$str" >> outfile_tsp_raw.log
				echo "$output" >> outfile_tsp_raw.log
				exit_status=$?
				if [[ $exit_status -eq 124 ]]; then
					echo TIMED OUT >> outfile_tsp.log
				elif [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then
					echo NO TIMING DATA >> outfile_tsp.log
				else
					echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile_tsp.log
				fi
			fi
		fi
	elif [[ $str == "ft70.4.sop" ]]; then
		echo "Starting..."
		start=1
	fi 
done

# join pairs of lines for processing by gnuplot
cat outfile_tsp.log | paste -d " " - - >> outfile_joined_tsp.log
