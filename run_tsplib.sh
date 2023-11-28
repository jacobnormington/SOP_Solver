#!/bin/bash

#if startInstance is set, it will add to the existing log files, picking up from a previous pause
startInstance=""					# the first instance to process, skipping any lexicographically earlier; if empty, begins at the first instance
endInstance=""						# the last instance to process, skipping any lexicographically later; if empty, begins at the first instance
dataset=0							# which instances to consider: 1 to exclude very hard instances, 0 for all instances, -1 for only very hard instances
num_threads=32 						# check 8, 16, 32



start=0 			# DO NOT USE, use startInstance instead
config="tsplib_config.txt"
if [[ $startInstance == "" ]]; then
	start=1
	rm -f outfile_tsp.log
	rm -f outfile_tsp_raw.log
	rm -f outfile_tsp_joined.log

	echo $num_threads threads >> outfile_tsp.log
	echo $num_threads threads >> outfile_tsp_raw.log
fi



# run the tests
for str in $(ls -v tsplib);
do
	if [[ $str == $startInstance ]]; then
		echo "Starting..."
		start=1
	fi

	if [[ $start == 1 ]]; then
		if [[ ${str: -4:4} == ".sop" ]]; then #discard anything other than sop instance files
			if [[ 	$dataset == 0 ||
					($dataset -gt "0" && ($str != "ESC78.sop" && $str != "ft53.1.sop" && $str != "ft53.2.sop" && $str != "ft53.3.sop" && $str != "ft70.1.sop" && $str != "ft70.2.sop" 
						&& $str != "ft70.3.sop" && $str != "kro124p.1.sop" && $str != "kro124p.2.sop" && $str != "kro124p.3.sop" && $str != "kro124p.4.sop" && $str != "p43.1.sop" 
						&& $str != "p43.2.sop" && $str != "prob.100.sop" && $str != "rbg323a.sop" && $str != "rbg341a.sop" && $str != "rbg358a.sop" && $str != "rbg378a.sop" 
						&& $str != "ry48p.1.sop" && $str != "ry48p.2.sop" && $str != "ry48p.3.sop")) ||
					($dataset -lt "0" && !($str != "ESC78.sop" && $str != "ft53.1.sop" && $str != "ft53.2.sop" && $str != "ft53.3.sop" && $str != "ft70.1.sop" && $str != "ft70.2.sop" 
						&& $str != "ft70.3.sop" && $str != "kro124p.1.sop" && $str != "kro124p.2.sop" && $str != "kro124p.3.sop" && $str != "kro124p.4.sop" && $str != "p43.1.sop" 
						&& $str != "p43.2.sop" && $str != "prob.100.sop" && $str != "rbg323a.sop" && $str != "rbg341a.sop" && $str != "rbg358a.sop" && $str != "rbg378a.sop" 
						&& $str != "ry48p.1.sop" && $str != "ry48p.2.sop" && $str != "ry48p.3.sop"))
				]]; then 
				echo $str >> outfile_tsp.log
				output=$(./sop_solver tsplib/$str $num_threads $config)
				echo "$str" >> outfile_tsp_raw.log
				echo "$output" >> outfile_tsp_raw.log

				if [[ $( echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" ) == "" ]]; then # time and cost not printed at the end
					echo NO TIMING DATA >> outfile_tsp.log
				else # output timing data to file
					echo $output | egrep -o "[[:digit:]]+,[[:digit:]]+(.[[:digit:]]+)?" >> outfile_tsp.log # add " | cut -d, -f2" to find only the time
				fi

				if [[ $( echo "$output" | egrep -o "instance timed out" ) != "" ]]; then # test times out
					echo TIMED OUT >> outfile_tsp.log
				fi

				# if [[ $( echo $output | egrep -o "Total Progress = 100%" ) == "" ]]; then # if test doesn't end in 100% progress
				# 	echo PROGRESS BAR NOT 100% >> outfile_tsp.log
				# 	echo PROGRESS BAR NOT 100% >> outfile_tsp_raw.log
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
# cat outfile_tsp.log | paste -d " " - - >> outfile_joined.log
